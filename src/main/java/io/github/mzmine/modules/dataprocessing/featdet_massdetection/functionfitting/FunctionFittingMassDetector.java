/*
 * Copyright (c) 2004-2022 The MZmine Development Team
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

package io.github.mzmine.modules.dataprocessing.featdet_massdetection.functionfitting;

import io.github.mzmine.datamodel.MassList;
import io.github.mzmine.datamodel.MassSpectrum;
import io.github.mzmine.datamodel.Scan;
import io.github.mzmine.modules.dataprocessing.featdet_massdetection.MassDetector;
import io.github.mzmine.parameters.ParameterSet;
import io.github.mzmine.util.collections.BinarySearch.DefaultTo;

import java.util.ArrayList;
import java.util.List;
import org.jetbrains.annotations.NotNull;
import org.apache.commons.math3.fitting.GaussianCurveFitter;
import org.apache.commons.math3.fitting.PolynomialCurveFitter;
import org.apache.commons.math3.fitting.WeightedObservedPoints;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

/**
 * This function searches for peaks in a mass spectrum using the Z-Score
 * algorithm.
 * 
 * Inspired from
 * https://stackoverflow.com/questions/22583391/peak-signal-detection-in-realtime-timeseries-data/22640362#22640362
 */
public class FunctionFittingMassDetector implements MassDetector {

  private final double threshold;
  private final int lag;
  private final int minPointsForPeak;
  private final boolean useMedian;
  private final double influence;
  private final double noiseLevel;

  public FunctionFittingMassDetector() {
    this(0, 0, 0, false, 0, 0);
  }

  public FunctionFittingMassDetector(
      final int lag,
      final double threshold,
      final int minPointsForPeak,
      final boolean useMedian,
      final double influence,
      final double noiseLevel) {

    this.lag = lag;
    this.threshold = threshold;
    this.minPointsForPeak = minPointsForPeak;
    this.useMedian = useMedian;
    this.influence = influence;
    this.noiseLevel = noiseLevel;
  }

  @Override
  public FunctionFittingMassDetector create(ParameterSet parameters) {
    // Read the paramters

    return new FunctionFittingMassDetector(
        parameters
            .getParameter(FunctionFittingMassDetectorParameters.lag)
            .getValue(),
        parameters
            .getParameter(FunctionFittingMassDetectorParameters.threshold)
            .getValue(),
        parameters
            .getParameter(FunctionFittingMassDetectorParameters.minPointsForPeak)
            .getValue(),
        parameters
            .getParameter(FunctionFittingMassDetectorParameters.useMedian)
            .getValue(),
        parameters
            .getParameter(FunctionFittingMassDetectorParameters.influence)
            .getValue(),
        parameters
            .getParameter(FunctionFittingMassDetectorParameters.noiseLevel)
            .getValue());
  }

  @Override
  public boolean filtersActive() {
    return true; // always active as z score only is used ?
  }

  @Override
  public double[][] getMassValues(MassSpectrum scan) {
    // Check the type of scan 
    // accept only the Scan 
    if (!(scan instanceof Scan)) {
      throw new IllegalArgumentException("The scan is not of type Scan");
    }

    return getMassValues((Scan) scan);
  }

  public double[][] getMassValues(Scan scan) {

    // Check that scan already has detected peaks
    MassList massList = scan.getMassList();
    if (massList == null) {
      throw new IllegalArgumentException("Mass list is required for function fitting. Please run another mass detector first.");
    }

    int numberOfMasses = massList.getNumberOfDataPoints();
    // Storage of the peaks
    List<Double> peaksMasses = new ArrayList<Double>();
    List<Double> peaksIntensities = new ArrayList<Double>();

    WeightedObservedPoints points = new WeightedObservedPoints();
    
    for (int i = 0; i < numberOfMasses; i++) {
      
    


      double mass = massList.getMzValue(i);
      // find the values in the scan around that mass 
      int thisPeakIndex = scan.binarySearch(mass, DefaultTo.CLOSEST_VALUE);

      int indexMin = thisPeakIndex - 2;
      int indexMax = thisPeakIndex + 2;

      if (indexMin < 0) {
        indexMin = 0;
      }
      if (indexMax >= scan.getNumberOfDataPoints()) {
        indexMax = scan.getNumberOfDataPoints() - 1;
      }

      if ((indexMax-indexMin) < 4) {
        // Cannot fit that peak
        // Simply put it in the output 
        peaksMasses.add(massList.getMzValue(i));
        peaksIntensities.add(massList.getIntensityValue(i));
        continue;
      }
    
      // add the points to the fitting
      for (int j = indexMin; j <= indexMax; j++) {
        points.add(scan.getMzValue(j), scan.getIntensityValue(j));
      }

      // init stats instance
        
      GaussianCurveFitter fitter = GaussianCurveFitter.create().withMaxIterations(1000);

      // Retrieve fitted parameters (coefficients of the polynomial function)
      try {
        final double[] coeff = fitter.fit(points.toList());
      
        double normalization = coeff[0];
        double mean = coeff[1];
        // double sigma = coeff[2];

        // Check that the value does not differ too much from the original value
        if (Math.abs(mean - mass) > 1.0) {
          // Could not fit that peak
          continue;
        }
  
        peaksMasses.add(mean);
        // Intensity is the value of the parabola at the apex
        peaksIntensities.add(normalization);
      } catch (Exception e) {
        // Could not fit that peak
        continue;
      }

    }


    // Assing the peaks to arrays of m/z and intensity for return
    int numberOfPeaks = peaksMasses.size();
    if (peaksIntensities.size() != numberOfPeaks) {
      // Should never happen
      throw new IllegalStateException("Number of peaks and intensities do not match");
    }
    double[] mzs = new double[numberOfPeaks];
    double[] intensities = new double[numberOfPeaks];

    for (int i = 0; i < numberOfPeaks; i++) {
      mzs[i] = peaksMasses.get(i);
      intensities[i] = peaksIntensities.get(i);
    }

    return new double[][] { mzs, intensities };
  }

  @Override
  public @NotNull String getName() {
    return "Function Fitting";
  }

  @Override
  public @NotNull Class<? extends ParameterSet> getParameterSetClass() {
    return FunctionFittingMassDetectorParameters.class;
  }

}
