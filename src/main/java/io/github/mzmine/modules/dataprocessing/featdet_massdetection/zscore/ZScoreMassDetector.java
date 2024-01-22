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

package io.github.mzmine.modules.dataprocessing.featdet_massdetection.zscore;

import io.github.mzmine.datamodel.MassSpectrum;
import io.github.mzmine.modules.dataprocessing.featdet_massdetection.MassDetector;
import io.github.mzmine.parameters.ParameterSet;
import java.util.ArrayList;
import java.util.List;
import org.jetbrains.annotations.NotNull;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

/**
 * This function searches for peaks in a mass spectrum using the Z-Score
 * algorithm.
 * 
 * Inspired from
 * https://stackoverflow.com/questions/22583391/peak-signal-detection-in-realtime-timeseries-data/22640362#22640362
 */
public class ZScoreMassDetector implements MassDetector {

  private final double threshold;
  private final int lag;
  private final int minPointsForPeak;
  private final boolean useMedian;
  private final double influence;
  private final double noiseLevel;

  public ZScoreMassDetector() {
    this(0, 0, 0, false, 0, 0);
  }

  public ZScoreMassDetector(
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
  public ZScoreMassDetector create(ParameterSet parameters) {
    // Read the paramters

    return new ZScoreMassDetector(
        parameters
            .getParameter(ZScoreMassDetectorParameters.lag)
            .getValue(),
        parameters
            .getParameter(ZScoreMassDetectorParameters.threshold)
            .getValue(),
        parameters
            .getParameter(ZScoreMassDetectorParameters.minPointsForPeak)
            .getValue(),
        parameters
            .getParameter(ZScoreMassDetectorParameters.useMedian)
            .getValue(),
        parameters
            .getParameter(ZScoreMassDetectorParameters.influence)
            .getValue(),
        parameters
            .getParameter(ZScoreMassDetectorParameters.noiseLevel)
            .getValue());
  }

  @Override
  public boolean filtersActive() {
    return true; // always active as z score only is used ?
  }

  @Override
  public double[][] getMassValues(MassSpectrum scan) {

    // Set up variables for the algorithm
    int scanSize = scan.getNumberOfDataPoints();

    // Make sure the lag value is okay
    if (lag < 2) {
      throw new IllegalArgumentException("The lag size: " + lag + " is smaller than 2");
    }
    if (scanSize < lag) {
      throw new IllegalArgumentException("The scan size: " + scanSize + " is smaller than the lag: " + lag);
    }

    if (minPointsForPeak < 1) {
      throw new IllegalArgumentException("The minPointsForPeak value: " + minPointsForPeak + " is smaller than 1");
    }

    // Storage of the peaks
    List<Integer> indexesPeaks = new ArrayList<Integer>();

    // init stats instance
    DescriptiveStatistics stats = new DescriptiveStatistics();

    // the average of the rolling lag
    double[] avgFilter = new double[scanSize];

    // the standard deviation of the rolling lag
    double[] stdFilter = new double[scanSize];
    // local zscore values
    double[] zscore = new double[scanSize];
    double[] filteredSignal = new double[scanSize];

    // Will say if the current position is peaking up
    int isPeaking = 0;

    // initialize the filter assuming no peak !
    for (int i = 0; i < lag; i++) {
      filteredSignal[i] = scan.getIntensityValue(i);
    }

    // init avgFilter and stdFilter
    for (int i = lag; i < scanSize; i++) {
      // Get the current intensity
      double thisIntensity = scan.getIntensityValue(i);

      if ((thisIntensity < noiseLevel) && (isPeaking == 0)) {
        // Noise signal, we don't care about it
        filteredSignal[i] = thisIntensity;
        continue;
      }
      // Create the local statisic
      // Todo: this is not efficient, we should update the stats instead of creating a
      // new one
      for (int j = 1; j <= lag; j++) {
        stats.addValue(filteredSignal[i - j]);
      }
      if (useMedian) {
        avgFilter[i] = stats.getPercentile(0.5);
      } else {
        avgFilter[i] = stats.getMean();
      }
      // getStandardDeviation() uses sample variance
      stdFilter[i] = Math.sqrt(stats.getPopulationVariance());
      
      // This is an adapted version of the z-score algorithm
      // It is adapted to work with a rolling mean and standard deviation
      zscore[i] = (thisIntensity - avgFilter[i]) / stdFilter[i];
      stats.clear();

      // Z-Score comparison: check if peaks stands out from the average
      if (zscore[i] > threshold) {
        // We are on a positive peaking zone
        isPeaking++;
        filteredSignal[i] = influence * thisIntensity + (1 - influence) * filteredSignal[i - 1];
      } else {
        // We are not on a positive peaking zone
        filteredSignal[i] = thisIntensity;
        // peaking stopped
        if (isPeaking > 0) {
          // Check the peak is wide enough
          if (isPeaking >= minPointsForPeak) {
            // We have a found a peaking zone, now we need to determine which position
            // within the zone
            // is the actual peak
            // We take the highest point in the zone
            double max = Double.MIN_VALUE;
            int maxIndex = -1;
            for (int j = i-isPeaking; j < i; j++) {
              double val = scan.getIntensityValue(j);
              if (val > max) {
                maxIndex = j;
                max = val;
              }
            }
            if (maxIndex == -1) {
              throw new IllegalStateException("maxIndex should not be -1");
            }
            indexesPeaks.add(maxIndex);
          }
          // Reset for the next interation
          isPeaking = 0;
        }
      }
    }

    // Assing the peaks to arrays of m/z and intensity for return
    int numberOfPeaks = indexesPeaks.size();
    double[] mzs = new double[numberOfPeaks];
    double[] intensities = new double[numberOfPeaks];

    int ind;
    for (int i = 0; i < numberOfPeaks; i++) {
      ind = indexesPeaks.get(i);
      mzs[i] = scan.getMzValue(ind);
      intensities[i] = scan.getIntensityValue(ind);
    }

    return new double[][] { mzs, intensities };
  }

  @Override
  public @NotNull String getName() {
    return "Z-Score";
  }

  @Override
  public @NotNull Class<? extends ParameterSet> getParameterSetClass() {
    return ZScoreMassDetectorParameters.class;
  }

}
