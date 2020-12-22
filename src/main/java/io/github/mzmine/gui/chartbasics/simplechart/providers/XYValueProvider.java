/*
 * Copyright 2006-2020 The MZmine Development Team
 *
 * This file is part of MZmine.
 *
 * MZmine is free software; you can redistribute it and/or modify it under the terms of the GNU
 * General Public License as published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 *
 * MZmine is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
 * the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
 * Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with MZmine; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
 * USA
 */

package io.github.mzmine.gui.chartbasics.simplechart.providers;

import io.github.mzmine.gui.chartbasics.simplechart.datasets.ColoredXYDataset;
import java.util.List;

/**
 * This interface is used to provide a dataset with x and y values. The amount of x and y values has
 * to be equal and is checked via the {@link List#size()} method.
 * <p></p>
 * The values are not grabbed during the creation of the dataset. After initialising the dataset
 * (e.g. {@link ColoredXYDataset}) a thread is started where the values of the dataset can be
 * calculated or loaded from disk. For that operation, the {@link XYValueProvider#computeValues()}
 * method is used. The implementing class can supply information on the progress of the operation
 * via the method {@link XYValueProvider#getComputationFinishedPercentage()}, which will be
 * represented in the task bar.
 * <p></p>
 * When the computation ({@link XYValueProvider#computeValues} has finished, the values are loaded
 * into the dataset via the {@link XYValueProvider#getDomainValues()} and {@link
 * XYValueProvider#getRangeValues()} methods.
 * <p></p>
 * After the dataset has been loaded successfully, the chart is automatically updated via a {@link
 * org.jfree.chart.JFreeChart#fireChartChanged()} event.
 *
 * @author https://github.com/SteffenHeu
 */
public interface XYValueProvider {

  /**
   * Called in a seperate thread to compute values or load them from disk after the dataset has been
   * created.
   */
  public void computeValues();

  /**
   * @return A sorted list of domain values. Index has to match the range value indices.
   */
  public List<Double> getDomainValues();

  /**
   * @return A sorted (ascending) list of range values. Index has to match the domain value indices.
   */
  public List<Double> getRangeValues();

  /**
   * Helper method to provide the user with progress information during {@link
   * XYValueProvider#computeValues()}.
   *
   * @return a finished percentage. (0.0-1.0)
   */
  public double getComputationFinishedPercentage();
}
