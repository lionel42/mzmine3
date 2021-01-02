/*
 *  Copyright 2006-2020 The MZmine Development Team
 *
 *  This file is part of MZmine.
 *
 *  MZmine is free software; you can redistribute it and/or modify it under the terms of the GNU
 *  General Public License as published by the Free Software Foundation; either version 2 of the
 *  License, or (at your option) any later version.
 *
 *  MZmine is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
 *  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
 *  Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License along with MZmine; if not,
 *  write to the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
 *  USA
 */

package io.github.mzmine.modules.visualization.rawdataoverviewims;

import com.google.common.collect.Range;
import io.github.mzmine.main.MZmineCore;
import io.github.mzmine.parameters.parametertypes.DoubleComponent;
import io.github.mzmine.parameters.parametertypes.ranges.DoubleRangeComponent;
import io.github.mzmine.parameters.parametertypes.selectors.ScanSelection;
import io.github.mzmine.parameters.parametertypes.selectors.ScanSelectionComponent;
import io.github.mzmine.parameters.parametertypes.tolerances.MZTolerance;
import io.github.mzmine.parameters.parametertypes.tolerances.MZToleranceComponent;
import java.text.NumberFormat;
import java.util.List;
import javafx.collections.FXCollections;
import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.control.Button;
import javafx.scene.control.Label;
import javafx.scene.control.ListCell;
import javafx.scene.control.ListView;
import javafx.scene.control.SelectionMode;
import javafx.scene.control.Tooltip;
import javafx.scene.layout.ColumnConstraints;
import javafx.scene.layout.FlowPane;
import javafx.scene.layout.GridPane;

public class IMSRawDataOverviewControlPanel extends GridPane {

  public static final String TOOLTIP_MOBILITYSCAN_NL = "Noise level for mobility scan processing. "
      + "Greatly impacts performance of this overview. Influences mobilogram, ion trace building "
      + "and the frame overview heatmap.";

  public static final String TOOLTIP_FRAME_NL = "Noise level for frame processing. Influences EIC"
      + " building and frame chart.";

  public static final String TOOLTIP_MZTOL = "m/z tolerance for EIC, ion trace and mobilogram "
      + "building";

  public static final String TOOLTIP_SCANSEL = "Scan selection for EIC, ion trace and mobilogram "
      + "building";


  private final IMSRawDataOverviewPane pane;
  private final NumberFormat mzFormat;
  private final NumberFormat intensityFormat;

  private MZTolerance mzTolerance;
  private ScanSelection scanSelection;
  private ListView<Range<Double>> mobilogramRangesList;

  private double frameNoiseLevel;
  private double mobilityScanNoiseLevel;

  IMSRawDataOverviewControlPanel(IMSRawDataOverviewPane pane, double frameNoiseLevel,
      double mobilityScanNoiseLevel, MZTolerance mzTolerance, ScanSelection scanSelection) {
    this.pane = pane;
    this.frameNoiseLevel = frameNoiseLevel;
    this.mobilityScanNoiseLevel = mobilityScanNoiseLevel;
    this.mzTolerance = mzTolerance;
    this.scanSelection = scanSelection;
    mzFormat = MZmineCore.getConfiguration().getMZFormat();
    intensityFormat = MZmineCore.getConfiguration().getIntensityFormat();
    initControlPanel();
  }

  private void initControlPanel() {
    DoubleComponent frameNoiseLevelComponent = new DoubleComponent(100, 0d, Double.MAX_VALUE,
        intensityFormat, frameNoiseLevel);
    frameNoiseLevelComponent.setText(intensityFormat.format(frameNoiseLevel));
    DoubleComponent mobilityScanNoiseLevelComponent = new DoubleComponent(100, 0d, Double.MAX_VALUE,
        intensityFormat, mobilityScanNoiseLevel);
    mobilityScanNoiseLevelComponent.setText(intensityFormat.format(mobilityScanNoiseLevel));
    MZToleranceComponent mzToleranceComponent = new MZToleranceComponent();
    mzToleranceComponent.setValue(mzTolerance);
    ScanSelectionComponent scanSelectionComponent = new ScanSelectionComponent();
    scanSelectionComponent.setValue(scanSelection);

    setPadding(new Insets(5));
    setVgap(5);
    getColumnConstraints().addAll(new ColumnConstraints(150),
        new ColumnConstraints());
    Label lblMobilityScanNoiseLevel = new Label("Mobility scan noise level");
    lblMobilityScanNoiseLevel.setTooltip(new Tooltip(TOOLTIP_MOBILITYSCAN_NL));
    add(lblMobilityScanNoiseLevel, 0, 0);
    add(mobilityScanNoiseLevelComponent, 1, 0);
    Label lblFrameNoiseLevel = new Label("Frame noise level");
    lblFrameNoiseLevel.setTooltip(new Tooltip(TOOLTIP_FRAME_NL));
    add(lblFrameNoiseLevel, 0, 1);
    add(frameNoiseLevelComponent, 1, 1);
    Label lblMzTol = new Label("m/z tolerance");
    lblMzTol.setTooltip(new Tooltip(TOOLTIP_MZTOL));
    add(lblMzTol, 0, 2);
    add(mzToleranceComponent, 1, 2);
    Label lblScanSel = new Label("Scan selection");
    lblScanSel.setTooltip(new Tooltip(TOOLTIP_SCANSEL));
    add(lblScanSel, 0, 3);
    add(scanSelectionComponent, 1, 3);

    DoubleRangeComponent mobilogramRangeComp = new DoubleRangeComponent(mzFormat);
    mobilogramRangesList = new ListView<>(
        FXCollections.observableArrayList());
    mobilogramRangesList.setMaxHeight(150);
    mobilogramRangesList.setMaxWidth(240);
    mobilogramRangesList.setPrefWidth(240);
    mobilogramRangesList.getSelectionModel().setSelectionMode(SelectionMode.SINGLE);
    Button addMzRange = new Button("Add range");
    addMzRange.setOnAction(e -> {
      Range<Double> range = mobilogramRangeComp.getValue();
      if (range == null) {
        return;
      }
      mobilogramRangesList.getItems().add(range);
    });

    Button removeMzRange = new Button("Remove range");
    removeMzRange.setOnAction(e -> mobilogramRangesList.getItems()
        .remove(mobilogramRangesList.getSelectionModel().getSelectedItem()));
    mobilogramRangesList.setCellFactory(param -> new ListCell<>() {
      @Override
      protected void updateItem(Range<Double> item, boolean empty) {
        super.updateItem(item, empty);
        if (empty || item == null) {
          setText(null);
          setGraphic(null);
          return;
        }
        setText(mzFormat.format(item.lowerEndpoint()) + " - " + mzFormat
            .format(item.upperEndpoint()));
        setGraphic(null);
      }
    });

    Button update = new Button("Update");
    update.setOnAction(e -> {
      try {
        frameNoiseLevel =
            Double.parseDouble(frameNoiseLevelComponent.getText());
        this.mobilityScanNoiseLevel =
            Double.parseDouble(mobilityScanNoiseLevelComponent.getText());
        frameNoiseLevelComponent.setText(intensityFormat.format(frameNoiseLevel));
        mobilityScanNoiseLevelComponent
            .setText(intensityFormat.format(this.mobilityScanNoiseLevel));
        scanSelection = scanSelectionComponent.getValue();
        mzTolerance = mzToleranceComponent.getValue();
        pane.setMzTolerance(mzTolerance);
        pane.setScanSelection(scanSelection);
        pane.setFrameNoiseLevel(frameNoiseLevel);
        pane.setMobilityScanNoiseLevel(this.mobilityScanNoiseLevel);

        pane.updateTicPlot();
        pane.onSelectedFrameChanged();
      } catch (NullPointerException | NumberFormatException ex) {
        ex.printStackTrace();
      }
    });

    add(new Label("EIC/Mobilogram ranges"), 0, 4);
    add(mobilogramRangesList, 1, 4, 1, 1);
    add(new Label("Range:"), 0, 5);
    add(mobilogramRangeComp, 1, 5);
    FlowPane buttons = new FlowPane(addMzRange, removeMzRange, update);
    buttons.setHgap(5);
    buttons.setAlignment(Pos.CENTER);
    add(buttons, 0, 6, 2, 1);
  }

  public MZTolerance getMzTolerance() {
    return mzTolerance;
  }

  public void setMzTolerance(MZTolerance mzTolerance) {
    this.mzTolerance = mzTolerance;
  }

  public ScanSelection getScanSelection() {
    return scanSelection;
  }

  public void setScanSelection(
      ScanSelection scanSelection) {
    this.scanSelection = scanSelection;
  }

  public List<Range<Double>> getMobilogramRangesList() {
    return mobilogramRangesList.getItems();
  }

  public double getFrameNoiseLevel() {
    return frameNoiseLevel;
  }

  public void setFrameNoiseLevel(double frameNoiseLevel) {
    this.frameNoiseLevel = frameNoiseLevel;
  }

  public double getMobilityScanNoiseLevel() {
    return mobilityScanNoiseLevel;
  }

  public void setMobilityScanNoiseLevel(double mobilityScanNoiseLevel) {
    this.mobilityScanNoiseLevel = mobilityScanNoiseLevel;
  }
}
