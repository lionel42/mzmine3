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

import java.text.NumberFormat;

import io.github.mzmine.main.MZmineCore;
import io.github.mzmine.modules.dataprocessing.featdet_massdetection.MassDetectorSetupDialog;
import io.github.mzmine.parameters.Parameter;
import io.github.mzmine.parameters.impl.SimpleParameterSet;
import io.github.mzmine.parameters.parametertypes.DoubleParameter;
import io.github.mzmine.parameters.parametertypes.IntegerParameter;
import io.github.mzmine.parameters.parametertypes.BooleanParameter;
import io.github.mzmine.util.ExitCode;

public class FunctionFittingMassDetectorParameters extends SimpleParameterSet {


  public static final IntegerParameter lag = new IntegerParameter(
      "Lag",
      "Size of the moving window that calculates the mean and standard deviation of historical data",
      7);
  public static final IntegerParameter minPointsForPeak = new IntegerParameter(
      "Points for peak",
      "Minimum number of points constituting a peak",
      2);

  public static final DoubleParameter threshold = new DoubleParameter(
      "Threshold",
      " The \"z-score\" at which the algorithm signals. ",
      NumberFormat.getNumberInstance(),
      3.5);

  public static final BooleanParameter useMedian = new BooleanParameter(
      "Use Median",
      "Use median instead of mean for calculating the z-score",
      false);

  public static final DoubleParameter influence = new DoubleParameter(
      "Influence",
      "The influence (between 0 and 1) of new signals on the mean and standard deviation. A higher influence means that signals will take longer to adapt to a new level.",
      NumberFormat.getNumberInstance(),
      0.5);

  public static final DoubleParameter noiseLevel = new DoubleParameter(
    "Noise level",
    "Intensities less than this value are interpreted as noise",
    MZmineCore.getConfiguration().getIntensityFormat(), 0.0);


  public FunctionFittingMassDetectorParameters() {
    super(new Parameter[] {lag, threshold, minPointsForPeak, useMedian, influence, noiseLevel},
        "https://mzmine.github.io/mzmine_documentation/module_docs/featdet_mass_detection/mass-detection-algorithms.html#wavelet-transform");
  }

  public ExitCode showSetupDialog(boolean valueCheckRequired) {
    MassDetectorSetupDialog dialog = new MassDetectorSetupDialog(
        valueCheckRequired, 
        FunctionFittingMassDetector.class,
        this);
    dialog.showAndWait();
    return dialog.getExitCode();
  }
}
