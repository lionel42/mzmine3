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

package io.github.mzmine.modules.dataanalysis.bubbleplots.logratioplot;

import io.github.mzmine.datamodel.data.Feature;
import io.github.mzmine.datamodel.data.FeatureList;
import io.github.mzmine.datamodel.data.FeatureListRow;
import java.util.Vector;
import java.util.logging.Logger;

import org.jfree.data.xy.AbstractXYZDataset;

import com.google.common.primitives.Doubles;

import io.github.mzmine.datamodel.RawDataFile;
import io.github.mzmine.modules.dataanalysis.bubbleplots.RTMZDataset;
import io.github.mzmine.parameters.ParameterSet;
import io.github.mzmine.util.MathUtils;
import io.github.mzmine.util.FeatureMeasurementType;

public class LogratioDataset extends AbstractXYZDataset implements RTMZDataset {

  /**
   * 
   */
  private static final long serialVersionUID = 1L;

  private Logger logger = Logger.getLogger(this.getClass().getName());

  private double[] xCoords = new double[0];
  private double[] yCoords = new double[0];
  private double[] colorCoords = new double[0];
  private FeatureListRow[] featureListRows = new FeatureListRow[0];

  private String datasetTitle;

  public LogratioDataset(FeatureList alignedFeatureList, ParameterSet parameters) {
    int numOfRows = alignedFeatureList.getNumberOfRows();

    RawDataFile groupOneFiles[] =
        parameters.getParameter(LogratioParameters.groupOneFiles).getValue();
    RawDataFile groupTwoFiles[] =
        parameters.getParameter(LogratioParameters.groupTwoFiles).getValue();
    FeatureMeasurementType measurementType =
        parameters.getParameter(LogratioParameters.measurementType).getValue();

    // Generate title for the dataset
    datasetTitle = "Logratio analysis";
    datasetTitle = datasetTitle.concat(" (");
    if (measurementType == FeatureMeasurementType.AREA)
      datasetTitle = datasetTitle.concat("Logratio of average feature areas");
    else
      datasetTitle = datasetTitle.concat("Logratio of average feature heights");
    datasetTitle = datasetTitle
        .concat(" in " + groupOneFiles.length + " vs. " + groupTwoFiles.length + " files");
    datasetTitle = datasetTitle.concat(")");
    logger.finest("Computing: " + datasetTitle);

    Vector<Double> xCoordsV = new Vector<Double>();
    Vector<Double> yCoordsV = new Vector<Double>();
    Vector<Double> colorCoordsV = new Vector<Double>();
    Vector<FeatureListRow> featureListRowsV = new Vector<FeatureListRow>();

    for (int rowIndex = 0; rowIndex < numOfRows; rowIndex++) {

      FeatureListRow row = alignedFeatureList.getRow(rowIndex);

      // Collect available feature intensities for selected files
      Vector<Double> groupOneFeatureIntensities = new Vector<Double>();
      for (int fileIndex = 0; fileIndex < groupOneFiles.length; fileIndex++) {
        Feature p = row.getFeature(groupOneFiles[fileIndex]);
        if (p != null) {
          if (measurementType == FeatureMeasurementType.AREA)
            groupOneFeatureIntensities.add((double) p.getArea());
          else
            groupOneFeatureIntensities.add((double) p.getHeight());
        }
      }
      Vector<Double> groupTwoFeatureIntensities = new Vector<Double>();
      for (int fileIndex = 0; fileIndex < groupTwoFiles.length; fileIndex++) {
        Feature p = row.getFeature(groupTwoFiles[fileIndex]);
        if (p != null) {
          if (measurementType == FeatureMeasurementType.AREA)
            groupTwoFeatureIntensities.add((double) p.getArea());
          else
            groupTwoFeatureIntensities.add((double) p.getHeight());
        }
      }

      // If there are at least one measurement from each group for this
      // feature then calc logratio and include this feature in the plot
      if ((groupOneFeatureIntensities.size() > 0) && (groupTwoFeatureIntensities.size() > 0)) {

        double[] groupOneInts = Doubles.toArray(groupOneFeatureIntensities);
        double groupOneAvg = MathUtils.calcAvg(groupOneInts);
        double[] groupTwoInts = Doubles.toArray(groupTwoFeatureIntensities);
        double groupTwoAvg = MathUtils.calcAvg(groupTwoInts);
        double logratio = Double.NaN;
        if (groupTwoAvg != 0.0)
          logratio = (double) (Math.log(groupOneAvg / groupTwoAvg) / Math.log(2.0));

        Double rt = (double) row.getAverageRT();
        Double mz = row.getAverageMZ();

        xCoordsV.add(rt);
        yCoordsV.add(mz);
        colorCoordsV.add(logratio);
        featureListRowsV.add(row);

      }

    }

    // Finally store all collected values in arrays
    xCoords = Doubles.toArray(xCoordsV);
    yCoords = Doubles.toArray(yCoordsV);
    colorCoords = Doubles.toArray(colorCoordsV);
    featureListRows = featureListRowsV.toArray(new FeatureListRow[0]);

  }

  public String toString() {
    return datasetTitle;
  }

  @Override
  public int getSeriesCount() {
    return 1;
  }

  @Override
  public Comparable<?> getSeriesKey(int series) {
    if (series == 0)
      return new Integer(1);
    else
      return null;
  }

  public Number getZ(int series, int item) {
    if (series != 0)
      return null;
    if ((colorCoords.length - 1) < item)
      return null;
    return colorCoords[item];
  }

  public int getItemCount(int series) {
    return xCoords.length;
  }

  public Number getX(int series, int item) {
    if (series != 0)
      return null;
    if ((xCoords.length - 1) < item)
      return null;
    return xCoords[item];
  }

  public Number getY(int series, int item) {
    if (series != 0)
      return null;
    if ((yCoords.length - 1) < item)
      return null;
    return yCoords[item];
  }

  public FeatureListRow getFeatureListRow(int item) {
    return featureListRows[item];
  }

}
