/*
 * Copyright (c) 2004-2024 The MZmine Development Team
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package io.github.mzmine.modules.io.import_rawdata_tw_hdf5;

import com.google.common.collect.Range;
import io.github.mzmine.datamodel.MZmineProject;
import io.github.mzmine.datamodel.MassSpectrumType;
import io.github.mzmine.datamodel.PolarityType;
import io.github.mzmine.datamodel.RawDataFile;
import io.github.mzmine.datamodel.Scan;
import io.github.mzmine.datamodel.features.SimpleFeatureListAppliedMethod;
import io.github.mzmine.datamodel.impl.SimpleScan;
import io.github.mzmine.modules.MZmineModule;
import io.github.mzmine.parameters.ParameterSet;
import io.github.mzmine.parameters.parametertypes.selectors.ScanSelection;
import io.github.mzmine.taskcontrol.AbstractTask;
import io.github.mzmine.taskcontrol.TaskStatus;
import io.github.mzmine.util.ExceptionUtils;
import io.github.mzmine.util.collections.BinarySearch;
import io.github.mzmine.util.collections.IndexRange;

import java.io.File;
import java.io.IOException;
import java.time.Instant;
import java.util.Hashtable;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.jetbrains.annotations.NotNull;
import ucar.ma2.Array;
import ucar.ma2.Index;
import ucar.nc2.Attribute;
import ucar.nc2.NetcdfFile;
import ucar.nc2.Variable;
import ucar.nc2.Group;

/**
 *
 */
public class TwHdf5ImportTask extends AbstractTask {

  private final Logger logger = Logger.getLogger(this.getClass().getName());

  private NetcdfFile inputFile;

  private int parsedScans;
  private int totalScans = 0, validScans = 0;

  private  int tofidStart = 0;
  private int tofidEnd;

  private Hashtable<Integer, Double> scansRetentionTimes;

  private final File file;
  private final MZmineProject project;
  private final RawDataFile newMZmineFile;
  private final ParameterSet parameters;
  private final Class<? extends MZmineModule> module;

  private Variable tofdataVariable, timingDataVariable;
  private Group fullSpectraGroup;

  private int massCalMode;
  private double[] massCalibrationParameters = new double[5];
  private double[] mzValues, mzValuesScans;

  private final int dimTime0 = 0;
  private final int dimTime1 = 1;
  // private final int dimSegment = 2;
  private final int dimTofid = 3;

  public TwHdf5ImportTask(MZmineProject project, File fileToOpen, RawDataFile newMZmineFile,
      @NotNull final Class<? extends MZmineModule> module, @NotNull final ParameterSet parameters,
      @NotNull Instant moduleCallDate) {
    super(null, moduleCallDate); // storage in raw data file
    this.project = project;
    this.file = fileToOpen;
    this.newMZmineFile = newMZmineFile;
    this.parameters = parameters;
    this.module = module;
  }

  /**
   * @see io.github.mzmine.taskcontrol.Task#getFinishedPercentage()
   */
  @Override
  public double getFinishedPercentage() {
    return totalScans == 0 ? 0 : (double) parsedScans / validScans;
  }

  /**
   * @see java.lang.Runnable#run()
   */
  @Override
  public void run() {

    // Update task status
    setStatus(TaskStatus.PROCESSING);
    logger.info("Started parsing file " + file);

    try {

      // Open file
      this.startReading();

      // Get the scan range to read
      ScanSelection scanSelection = this.parameters.getParameter(TwHdf5ImportParameters.scanSelection).getValue();
      Range<Integer> scanRangeInt = scanSelection.getScanNumberRange();
      Range<Double> retentionTimeRange = scanSelection.getScanRTRange();

      // Parse scans
      for (int i = 0; i < scansRetentionTimes.size(); i++) {
        if (isCanceled()) {
          return;
        }
        Double retentionTime = scansRetentionTimes.get(i);
        if (retentionTime == null) {
          throw new IOException("Could not read retention time for scan " + i);
        }
        // Add the scan only if it is in range
        if ((scanRangeInt == null && retentionTimeRange == null)
            || (scanRangeInt != null
                && scanRangeInt.contains(i))
            || (retentionTimeRange != null
                && retentionTimeRange.contains(retentionTime))) {
          newMZmineFile.addScan(this.readScan(i));
        }
        parsedScans++;

      }

      // Close file
      this.finishReading();
      newMZmineFile.getAppliedMethods()
          .add(new SimpleFeatureListAppliedMethod(module, parameters, getModuleCallDate()));
      project.addFile(newMZmineFile);

    } catch (Throwable e) {
      logger.log(Level.SEVERE, "Could not open file " + file.getPath(), e);
      setErrorMessage(ExceptionUtils.exceptionToString(e));
      setStatus(TaskStatus.ERROR);
      return;
    }

    logger.info("Finished parsing " + file + ", parsed " + parsedScans + " scans");

    // Update task status
    setStatus(TaskStatus.FINISHED);

  }

  @Override
  public String getTaskDescription() {
    return "Opening file " + file;
  }

  public void startReading() throws IOException {

    // Open NetCDF-file
    try {
      inputFile = NetcdfFile.open(file.getPath());
    } catch (Exception e) {
      logger.severe(e.toString());
      throw (new IOException("Couldn't open input file" + file));
    }

    /*
     * DEBUG: dump all variables for (Variable v : inputFile.getVariables()) {
     * System.out.println("variable " + v.getShortName()); }
     */

    // Find mass_values and intensity_values variables
    tofdataVariable = inputFile.findVariable("FullSpectra/TofData");
    if (tofdataVariable == null) {
      logger.severe("Could not find variable FullSpectra/TofDatas");
      throw (new IOException("Could not find variable mass_values"));
    }
    assert (tofdataVariable.getRank() == 4);

    // Read the default mass calibration parameters
    fullSpectraGroup = inputFile.findGroup("FullSpectra");

    // Get the mode of the default mass calibration
    Attribute attribute = fullSpectraGroup.findAttribute("MassCalibMode");
    if (attribute == null) {
      logger.severe("Could not find attribute MassCalibMode");
      throw (new IOException("Could not find attribute MassCalibMode"));
    }
    massCalMode = attribute.getNumericValue().intValue();

    // Loop over the mass calibration parameters to read them
    for (int i = 0; i < 5; i++) {
      // mode 0 and 1 have only 2 params
      if (massCalMode < 2 && i > 1) {
        break;
      }
      // mode 2 has only 3 params
      if (massCalMode == 2 && i > 2) {
        break;
      }
      // mode 3 has only 4 params
      if (massCalMode == 3 && i > 3) {
        break;
      }

      // Create the attribute name
      String attributeName = "MassCalibration_p" + (i + 1);
      fullSpectraGroup = inputFile.findGroup("FullSpectra");
      Attribute attributeParam = fullSpectraGroup.findAttribute(attributeName);
      if (attributeParam == null) {
        logger.severe("Could not find attribute " + attributeName);
        throw (new IOException("Could not find attribute " + attributeName));
      }
      massCalibrationParameters[i] = attributeParam.getNumericValue().doubleValue();
    }

    // Generate the mass axis
    mzValues = new double[tofdataVariable.getDimension(dimTofid).getLength()];
    switch (massCalMode) {
      case 0:

        for (int i = 0; i < mzValues.length; i++) {
          double intermediate = ((i - massCalibrationParameters[1]) / massCalibrationParameters[0]);
          mzValues[i] = intermediate * intermediate;
        }
        break;

      case 1:
        for (int i = 0; i < mzValues.length; i++) {
          double intermediate = (massCalibrationParameters[0] / (i - massCalibrationParameters[1]));
          mzValues[i] = intermediate * intermediate;
        }
        break;
      case 2:
        double power = (1 / massCalibrationParameters[2]);
        for (int i = 0; i < mzValues.length; i++) {
          double intermediate = ((i - massCalibrationParameters[1]) / massCalibrationParameters[0]);
          mzValues[i] = Math.pow(intermediate, power);
        }
        break;

      default:
        throw new IllegalArgumentException("Unexpected/Unimplemented value for mass calibration mode: " + massCalMode);
    }

    // Find the mz range specified  
    Range<Double> mzRange = parameters.getParameter(TwHdf5ImportParameters.mzRange).getValue();
    
    if (mzRange != null) {
      // search where the mz value should go in the mzValues array
      IndexRange tofidRange = BinarySearch.indexRange(mzValues, mzRange);
      tofidStart = tofidRange.min();
      tofidEnd = tofidRange.maxExclusive();
      // Crop the mzValues array
      mzValuesScans = new double[tofidEnd - tofidStart];
      System.arraycopy(mzValues, tofidStart, mzValuesScans, 0, mzValuesScans.length);
     } else {
      tofidStart = 0;
      tofidEnd = mzValues.length;
      mzValuesScans = mzValues;
    }

    timingDataVariable = inputFile.findVariable("TimingData/BufTimes");
    if (timingDataVariable == null) {
      logger.severe("Could not find variable TimingData/BufTimes");
      throw (new IOException("Could not find variable TimingData/BufTimes"));
    }
    assert (timingDataVariable.getRank() == 2);

    // Read the timing data values (over the 2 timing dimensions)
    totalScans = timingDataVariable.getDimension(dimTime0).getLength() * timingDataVariable
        .getDimension(dimTime1).getLength();

    

    Array scanTimeArray = null;
    try {
      scanTimeArray = timingDataVariable.read();
    } catch (Exception e) {
      logger.severe(e.toString());
      throw (new IOException(
          "Could not read from variable TimingData/BufTimes from file " + file));
    }

    // Check the size of the array
    if (scanTimeArray.getSize() != totalScans) {
      logger.severe("Size of array TimingData/BufTimes is not equal to the number of scans");
      throw (new IOException(
          "Size of array TimingData/BufTimes is not equal to the number of scans"));
    }

    // Collect information about retention times, start positions and
    // lengths for scans
    scansRetentionTimes = new Hashtable<Integer, Double>();
    double retentionTime = -1;
    for (int i = 0; i < totalScans; i++) {

      Integer scanNum = i;

      // int dim0 = i / timingDataVariable.getDimension(dimTime1).getLength();
      // int dim1 = i % timingDataVariable.getDimension(dimTime1).getLength();
      double newRetentionTime = scanTimeArray.getDouble(i);
      if ((newRetentionTime < retentionTime) || (newRetentionTime - retentionTime > 100000)) {
        // This wierd things can happen at the end of the file
        break;
      }
      validScans++;
      scansRetentionTimes.put(scanNum, newRetentionTime);
      retentionTime = newRetentionTime;

    }

  }

  public void finishReading() throws IOException {
    inputFile.close();
  }

  /**
   * Reads one scan from the file. Requires that general information has already
   * been read.
   */
  private Scan readScan(int scanNum) throws IOException {

    // End of file
    if (scanNum >= validScans) {
      throw (new IOException("End of file reached"));
    }

    // Get retention time of the scan
    Float retentionTime = scansRetentionTimes.get(scanNum).floatValue();

    // What are these variables for?
    PolarityType polarity = PolarityType.UNKNOWN;
    String scanDefinition = "";

    // Get coordinates in the file
    int dim0 = scanNum / timingDataVariable.getDimension(dimTime1).getLength();
    int dim1 = scanNum % timingDataVariable.getDimension(dimTime1).getLength();

    int[] origin = new int[] { dim0, dim1, 0, tofidStart };
    int[] size = new int[] { 1, 1, 1, tofidEnd - tofidStart };

    // Read mass and intensity values
    Array intensityValueArray;
    try {
      intensityValueArray = tofdataVariable.read(origin, size);
    } catch (Exception e) {
      logger.log(Level.SEVERE, "Could not read from tofdataVariable.",
          e);
      throw (new IOException("Could not read from tofdataVariable."));
    }

    // prepare variables for extraction
    Index intensityValuesIndex = intensityValueArray.getIndex();
    int arrayLength = intensityValueArray.getShape()[dimTofid];
    double[] intensityValues = new double[arrayLength];

    for (int j = 0; j < arrayLength; j++) {
      intensityValues[j] = intensityValueArray.getDouble(intensityValuesIndex.set(0, 0, 0, j));
    }

    scanNum++;

    SimpleScan buildingScan = new SimpleScan(
        newMZmineFile, scanNum, 1, retentionTime, null,
        mzValuesScans, intensityValues,
        MassSpectrumType.PROFILE, polarity, scanDefinition, null);

    return buildingScan;

  }

}
