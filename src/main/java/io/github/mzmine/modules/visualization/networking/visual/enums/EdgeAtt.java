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

package io.github.mzmine.modules.visualization.networking.visual.enums;

public enum EdgeAtt implements GraphElementAttr {

  NONE, SCORE, SIM_N, TYPE, LABEL, GNPS_SCORE, DIFF_SCORE, SIM_SCORE, DIFF_N, DELTA_MZ, NUMBER_OF_COLLAPSED_EDGES, NEIGHBOR_DISTANCE;

  @Override
  public String toString() {
    return super.toString().replaceAll("_", " ");
  }

  public boolean isNumber() {
    return switch (this) {
      case TYPE, LABEL, NONE -> false;
      case SCORE, SIM_N, GNPS_SCORE, DIFF_SCORE, SIM_SCORE, DIFF_N, DELTA_MZ, NUMBER_OF_COLLAPSED_EDGES, NEIGHBOR_DISTANCE ->
          true;
    };
  }

  @Override
  public boolean isReversed() {
    return this==NEIGHBOR_DISTANCE;
  }

  @Override
  public boolean isChangingDynamically() {
    return NEIGHBOR_DISTANCE == this;
  }

  @Override
  public GraphObject getGraphObject() {
    return GraphObject.EDGE;
  }
}
