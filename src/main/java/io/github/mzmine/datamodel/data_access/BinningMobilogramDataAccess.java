/*
 * Copyright 2006-2021 The MZmine Development Team
 *
 * This file is part of MZmine.
 *
 * MZmine is free software; you can redistribute it and/or modify it under the terms of the GNU
 * General Public License as published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 *
 * MZmine is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
 * the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with MZmine; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 */

package io.github.mzmine.datamodel.data_access;

import com.google.common.collect.Range;
import com.google.common.primitives.Booleans;
import io.github.mzmine.datamodel.Frame;
import io.github.mzmine.datamodel.IMSRawDataFile;
import io.github.mzmine.datamodel.MobilityScan;
import io.github.mzmine.datamodel.MobilityType;
import io.github.mzmine.datamodel.featuredata.IntensitySeries;
import io.github.mzmine.datamodel.featuredata.IonMobilitySeries;
import io.github.mzmine.datamodel.featuredata.MobilitySeries;
import io.github.mzmine.datamodel.featuredata.impl.SummedIntensityMobilitySeries;
import io.github.mzmine.datamodel.features.FeatureList.FeatureListAppliedMethod;
import io.github.mzmine.datamodel.features.ModularFeatureList;
import io.github.mzmine.main.MZmineCore;
import io.github.mzmine.modules.dataprocessing.featdet_ionmobilitytracebuilder.AdvancedImsTraceBuilderParameters;
import io.github.mzmine.modules.dataprocessing.featdet_ionmobilitytracebuilder.IonMobilityTraceBuilderModule;
import io.github.mzmine.modules.dataprocessing.featdet_ionmobilitytracebuilder.IonMobilityTraceBuilderParameters;
import io.github.mzmine.modules.dataprocessing.featdet_mobilogram_summing.MobilogramBinningModule;
import io.github.mzmine.modules.dataprocessing.featdet_mobilogram_summing.MobilogramBinningParameters;
import io.github.mzmine.modules.dataprocessing.featdet_recursiveimsbuilder.RecursiveIMSBuilderAdvancedParameters;
import io.github.mzmine.modules.dataprocessing.featdet_recursiveimsbuilder.RecursiveIMSBuilderModule;
import io.github.mzmine.modules.dataprocessing.featdet_recursiveimsbuilder.RecursiveIMSBuilderParameters;
import io.github.mzmine.parameters.ParameterSet;
import io.github.mzmine.util.IonMobilityUtils;
import io.github.mzmine.util.MemoryMapStorage;
import java.nio.DoubleBuffer;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.logging.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

/**
 * Used to efficiently access mobilogram data of a raw data file. The data can be binned by mobility
 * to generate less noisy mobilograms.
 *
 * @author https://github.com/SteffenHeu
 */
public class BinningMobilogramDataAccess implements IntensitySeries, MobilitySeries {

  private static final double MOBILITY_EPSILON = 0.00001;

  private static Logger logger = Logger.getLogger(BinningMobilogramDataAccess.class.getName());

  private final IMSRawDataFile dataFile;

  private final double[] intensities;
  private final double[] tempMobilities;
  private final double[] tempIntensities;
  private final double[] mobilities;
  private final double[] upperBinLimits;
  private final int binWidth;

  private final double approximateBinSize;

  public BinningMobilogramDataAccess(@NotNull final IMSRawDataFile rawDataFile,
      final int binWidth) {
    if (binWidth < 1) {
      throw new IllegalArgumentException("Illegal bin width (" + binWidth + ")");
    }
    dataFile = rawDataFile;

    final Map<Frame, Range<Double>> ranges = IonMobilityUtils.getUniqueMobilityRanges(rawDataFile);
    // multiple mobility ranges are possible in tims
    final int maxMobilityScans =
        rawDataFile.getFrames().stream().mapToInt(Frame::getNumberOfMobilityScans).max()
            .orElseThrow() * ranges.size();

    tempIntensities = new double[maxMobilityScans];
    tempMobilities = new double[maxMobilityScans];
    this.binWidth = binWidth;

    var entries = ranges.entrySet().stream().toList();
    List<Double> distinctMobilities = new ArrayList<>();
    final MobilityType mt = rawDataFile.getMobilityType();

    // find all possible mobility values
    for (int j = 0; j < entries.size(); j++) {
      final Entry<Frame, Range<Double>> entry = entries.get(j);
      final Frame frame = entry.getKey();

      for (final MobilityScan scan : frame.getMobilityScans()) {
        if (!distinctMobilities.isEmpty()
            // either not tims and current mobility > highest mobility
            && ((distinctMobilities.get(distinctMobilities.size() - 1) > scan.getMobility()
            && mt != MobilityType.TIMS)
            // or tims and current mobility < lowest mobility
            || (distinctMobilities.get(distinctMobilities.size() - 1) < scan.getMobility()
            && mt == MobilityType.TIMS))) {
          continue;
        }
        distinctMobilities.add(scan.getMobility());
      }
    }

    distinctMobilities.sort(Double::compare);

    final List<Double> upperLimits = new ArrayList<>();
    final List<Double> centerBins = new ArrayList<>();
    for (int i = 0; i < distinctMobilities.size(); i += binWidth) {

      int currentBins = 0;
      double summedMobility = 0d;
      for (int j = 0; i + j < distinctMobilities.size() && j < binWidth; j++) {
        summedMobility += distinctMobilities.get(i + j);
        currentBins++;
        double d = distinctMobilities.get(i + j);
//        logger.finest(() -> "Adding " + d + " to bin");
      }
      centerBins.add(summedMobility / Math.max(currentBins, 1));
      upperLimits.add(
          distinctMobilities.get(Math.min(i + binWidth - 1, distinctMobilities.size() - 1))
              + MOBILITY_EPSILON);
    }

    mobilities = centerBins.stream().mapToDouble(Double::doubleValue).toArray();
    upperBinLimits = upperLimits.stream().mapToDouble(Double::doubleValue).toArray();
    intensities = new double[mobilities.length];

    double previous = mobilities[0];
    double deltas = 0;
    for (int i = 1; i < mobilities.length; i++) {
      deltas += mobilities[i] - previous;
      previous = mobilities[i];
    }

    approximateBinSize = deltas / (mobilities.length - 2);
    logger.finest(
        () -> "Bin width set to " + binWidth + " scans. (approximately " + approximateBinSize + " "
            + rawDataFile.getMobilityType().getUnit() + ")");
  }

  @NotNull
  public static int getRecommendedBinWidth(IMSRawDataFile file) {
    final Frame frame = file.getFrame(0);
    switch (frame.getMobilityType()) {
      case NONE, DRIFT_TUBE, TRAVELING_WAVE, FAIMS -> {
        return 1;
      }
      case TIMS -> {
        int index = frame.getNumberOfMobilityScans() / 2;
        final double mob1 = frame.getMobilityScan(index).getMobility();
        final double mob2 = frame.getMobilityScan(index + 1).getMobility();
        final double delta = Math.abs(mob1 - mob2);
        return (int) Math.max(1d, 0.0008 / delta);
      }
    }
    return 1;
  }

  @Nullable
  public static Integer getPreviousBinningWith(@NotNull final ModularFeatureList flist,
      MobilityType mt) {
    List<FeatureListAppliedMethod> methods = flist.getAppliedMethods();

    Integer binWidth = null;
    for (int i = methods.size() - 1; i >= 0; i--) {
      FeatureListAppliedMethod method = methods.get(i);
      if (method.getModule()
          .equals(MZmineCore.getModuleInstance(IonMobilityTraceBuilderModule.class))) {
        final ParameterSet parameterSet = method.getParameters();
        final var advancedParam = parameterSet
            .getParameter(IonMobilityTraceBuilderParameters.advancedParameters).getValue();
        binWidth = switch (mt) {
          case TIMS ->
              advancedParam.getParameter(AdvancedImsTraceBuilderParameters.timsBinningWidth)
                  .getValue() ? advancedParam
                  .getParameter(AdvancedImsTraceBuilderParameters.timsBinningWidth)
                  .getEmbeddedParameter().getValue()
                  : MobilogramBinningParameters.DEFAULT_TIMS_BIN_WIDTH;
          case DRIFT_TUBE ->
              advancedParam.getParameter(AdvancedImsTraceBuilderParameters.dtimsBinningWidth)
                  .getValue() ? advancedParam
                  .getParameter(AdvancedImsTraceBuilderParameters.dtimsBinningWidth)
                  .getEmbeddedParameter().getValue()
                  : MobilogramBinningParameters.DEFAULT_DTIMS_BIN_WIDTH;
          case TRAVELING_WAVE ->
              advancedParam.getParameter(AdvancedImsTraceBuilderParameters.twimsBinningWidth)
                  .getValue() ? advancedParam
                  .getParameter(AdvancedImsTraceBuilderParameters.twimsBinningWidth)
                  .getEmbeddedParameter().getValue()
                  : MobilogramBinningParameters.DEFAULT_TWIMS_BIN_WIDTH;
          default -> null;
        };
        break;
      }

      if (method.getModule()
          .equals(MZmineCore.getModuleInstance(RecursiveIMSBuilderModule.class))) {
        final ParameterSet parameterSet = method.getParameters();
        final var advancedParam = parameterSet
            .getParameter(RecursiveIMSBuilderParameters.advancedParameters).getValue();
        binWidth = switch (mt) {
          case TIMS ->
              advancedParam.getParameter(RecursiveIMSBuilderAdvancedParameters.timsBinningWidth)
                  .getValue() ? advancedParam
                  .getParameter(RecursiveIMSBuilderAdvancedParameters.timsBinningWidth)
                  .getEmbeddedParameter().getValue()
                  : MobilogramBinningParameters.DEFAULT_TIMS_BIN_WIDTH;
          case DRIFT_TUBE ->
              advancedParam.getParameter(RecursiveIMSBuilderAdvancedParameters.dtimsBinningWidth)
                  .getValue() ? advancedParam
                  .getParameter(RecursiveIMSBuilderAdvancedParameters.dtimsBinningWidth)
                  .getEmbeddedParameter().getValue()
                  : MobilogramBinningParameters.DEFAULT_DTIMS_BIN_WIDTH;
          case TRAVELING_WAVE ->
              advancedParam.getParameter(RecursiveIMSBuilderAdvancedParameters.twimsBinningWidth)
                  .getValue() ? advancedParam
                  .getParameter(RecursiveIMSBuilderAdvancedParameters.twimsBinningWidth)
                  .getEmbeddedParameter().getValue()
                  : MobilogramBinningParameters.DEFAULT_TWIMS_BIN_WIDTH;
          default -> null;
        };
        break;
      }

      if (method.getModule().equals(MZmineCore.getModuleInstance(MobilogramBinningModule.class))) {
        final ParameterSet parameterSet = method.getParameters();
        binWidth = switch (mt) {
          case TIMS -> parameterSet.getParameter(MobilogramBinningParameters.timsBinningWidth)
              .getValue();
          case DRIFT_TUBE -> parameterSet
              .getParameter(MobilogramBinningParameters.dtimsBinningWidth).getValue();
          case TRAVELING_WAVE -> parameterSet
              .getParameter(MobilogramBinningParameters.twimsBinningWidth).getValue();
          default -> null;
        };
        break;
      }
    }
    if (binWidth == null) {
      logger.info(
          () -> "Previous binning width not recognised. Has the mobility type been implemented?");
      binWidth = getRecommendedBinWidth((IMSRawDataFile) flist.getRawDataFile(0));
    }
    return binWidth;
  }

  private void clearIntensities() {
    Arrays.fill(intensities, 0d);
  }

  /**
   * Re-bins an already summed mobilogram. Note that re-binning an already binned mobilogram with a
   * lower binnign width than before will lead to 0-intensity values. Consider using {@link
   * #setMobilogram(List)} instead.
   *
   * @param summedMobilogram
   */
  public void setMobilogram(@NotNull final SummedIntensityMobilitySeries summedMobilogram) {
    clearIntensities();

    final int numValues = summedMobilogram.getNumberOfValues();
    assert numValues <= tempIntensities.length;

    summedMobilogram.getIntensityValues(tempIntensities);
    summedMobilogram.getMobilityValues(tempMobilities);

    int rawIndex = 0;
    for (int binnedIndex = 0; binnedIndex < intensities.length && rawIndex < numValues;
        binnedIndex++) {
      // ensure we are above the current lower-binning-limit
      while (rawIndex < numValues && Double.compare(tempMobilities[rawIndex],
          binnedIndex == 0 ? 0d : upperBinLimits[binnedIndex - 1]) == -1) {
        rawIndex++;
      }

      // ensure we are below the current upper-binning-limit
      while (rawIndex < numValues
          && Double.compare(tempMobilities[rawIndex], upperBinLimits[binnedIndex]) == -1) {
        intensities[binnedIndex] += tempIntensities[rawIndex];
        rawIndex++;
      }
    }
  }

  /**
   * Constructs a binned summed mobilogram from the supplied list of individual mobilograms.
   *
   * @param mobilograms The list of {@link IonMobilitySeries}.
   */
  public void setMobilogram(@NotNull final List<IonMobilitySeries> mobilograms) {
    clearIntensities();

    int order = 1;
    if (!mobilograms.isEmpty()) {
      order =
          mobilograms.get(0).getSpectrum(0).getFrame().getMobilityType() == MobilityType.TIMS ? -1
              : +1;
    }

    for (IonMobilitySeries ims : mobilograms) {
      final int numValues = ims.getNumberOfValues();
      ims.getIntensityValues(tempIntensities);

      for (int i = 0; i < numValues; i++) {
        tempMobilities[i] = ims.getMobility(i);
      }

      // in tims, the mobilograms are sorted by decreasing order
      final int start = order == 1 ? 0 : numValues - 1;
      int rawIndex = start;

      boolean[] assigned = new boolean[numValues];
      Arrays.fill(assigned, false);

      for (int i = 0; i < upperBinLimits.length && rawIndex >= 0; i++) {

        // waters records DT = 0, so it cannot be 0
        final double binStart = i == 0 ? -MOBILITY_EPSILON : upperBinLimits[i - 1];
        final double binEnd = upperBinLimits[i];
        // if we are in the correct bin, add all values that fit
        while (rawIndex >= 0 && rawIndex < numValues && tempMobilities[rawIndex] < binEnd
            && tempMobilities[rawIndex] > binStart) {

          intensities[i] += tempIntensities[rawIndex];
          assigned[rawIndex] = true;
          rawIndex += order;
        }
      }

      long numAssigned = Booleans.asList(assigned).stream().filter(val -> val == true).count();
      if (numAssigned != numValues) {
        logger.info("assiged " + numAssigned + "/" + numValues);
      }
    }
  }

  public SummedIntensityMobilitySeries toSummedMobilogram(@Nullable MemoryMapStorage storage) {
    int firstNonZero = -1;
    int lastNonZero = -1;

    for (int i = 0; i < mobilities.length; i++) {
      if (firstNonZero == -1 && intensities[i] > 0d) {
        firstNonZero = Math.max(0, i - 1);
      }
      if (intensities[i] > 0d) {
        lastNonZero = i;
      }
    }
    firstNonZero = Math.max(firstNonZero, 0);
    lastNonZero = Math.min(lastNonZero + 1, mobilities.length - 1);

    return new SummedIntensityMobilitySeries(storage,
        Arrays.copyOfRange(mobilities, firstNonZero, lastNonZero),
        Arrays.copyOfRange(intensities, firstNonZero, lastNonZero));
  }

  @Override
  public DoubleBuffer getIntensityValueBuffer() {
    throw new UnsupportedOperationException(
        "This data access is designed to loop over intensities/mobilities.");
  }

  /**
   * @param dst a buffer to copy the intensities to. must be of appropriate size. If null is passed,
   *            the intensity array is returned directly. Do not modify.
   * @return The intensity values.
   */
  @Override
  public double[] getIntensityValues(double[] dst) {
    if (dst != null) {
      assert dst.length >= getNumberOfValues();
      System.arraycopy(intensities, 0, dst, 0, getNumberOfValues());
      return dst;
    } else {
      return intensities;
    }
  }

  public double[] getIntensityValues() {
    return intensities;
  }

  public double[] getMobilityValues() {
    return mobilities;
  }

  @Override
  public double getIntensity(int index) {
    return intensities[index];
  }

  @Override
  public int getNumberOfValues() {
    return intensities.length;
  }

  @Override
  public double getMobility(int index) {
    return mobilities[index];
  }

  /**
   * @param dst a buffer to copy the mobilities to. must be of appropriate size. If null is passed,
   *            the intensity array is returned directly. Do not modify.
   * @return The intensity values.
   */
  public double[] getMobilityValues(@Nullable double[] dst) {
    if (dst != null) {
      assert dst.length >= mobilities.length;
      System.arraycopy(mobilities, 0, dst, 0, mobilities.length);
      return dst;
    } else {
      return mobilities;
    }
  }

  public IMSRawDataFile getDataFile() {
    return dataFile;
  }

  public double getBinWidth() {
    return (double) binWidth;
  }

  /**
   * @return The approximate bin size in mobility units with respect to the raw data file.
   */
  public double getApproximateBinSize() {
    return approximateBinSize;
  }

}
