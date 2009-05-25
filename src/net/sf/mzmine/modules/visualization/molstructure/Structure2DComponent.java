/*
 * Copyright 2006-2009 The MZmine 2 Development Team
 * 
 * This file is part of MZmine 2.
 * 
 * MZmine 2 is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 2 of the License, or (at your option) any later
 * version.
 * 
 * MZmine 2 is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along with
 * MZmine 2; if not, write to the Free Software Foundation, Inc., 51 Franklin St,
 * Fifth Floor, Boston, MA 02110-1301 USA
 */

package net.sf.mzmine.modules.visualization.molstructure;

import java.awt.Point;
import java.io.StringReader;
import java.text.NumberFormat;
import java.util.Iterator;

import javax.swing.JLabel;
import javax.swing.JViewport;
import javax.swing.event.ChangeEvent;

import net.sf.mzmine.main.MZmineCore;

import org.openscience.cdk.ChemModel;
import org.openscience.cdk.applications.jchempaint.JChemPaintEditorPanel;
import org.openscience.cdk.applications.jchempaint.JChemPaintModel;
import org.openscience.cdk.controller.Controller2DModel;
import org.openscience.cdk.controller.PopupController2D;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemModel;
import org.openscience.cdk.io.MDLV2000Reader;
import org.openscience.cdk.renderer.Renderer2DModel;
import org.openscience.cdk.tools.MFAnalyser;
import org.openscience.cdk.tools.manipulator.ChemModelManipulator;

public class Structure2DComponent extends JChemPaintEditorPanel {

	public static final NumberFormat massFormater = MZmineCore.getMZFormat();
	private JLabel statusLabel;

	public Structure2DComponent(String structure) throws CDKException {
		this(structure, null);
	}

	public Structure2DComponent(String structure, JLabel statusLabel)
			throws CDKException {

		this.statusLabel = statusLabel;

		setShowStatusBar(false);
		setShowInsertTextField(false);
		setShowMenuBar(false);
		setShowToolBar(false);
		setOpaque(false);

		StringReader reader = new StringReader(structure);
		MDLV2000Reader molReader = new MDLV2000Reader(reader);
		ChemModel chemModel = new ChemModel();
		chemModel = (ChemModel) molReader.read(chemModel);
		processChemModel(chemModel);

		JChemPaintModel jcpModel = new JChemPaintModel();
		registerModel(jcpModel);
		setJChemPaintModel(jcpModel, null);

		Controller2DModel controller = jcpModel.getControllerModel();
		controller.setAutoUpdateImplicitHydrogens(true);
		controller.setDrawMode(Controller2DModel.LASSO);
		controller.setMovingAllowed(false);
		controller.setAutoUpdateImplicitHydrogens(true);

		Renderer2DModel renderer = jcpModel.getRendererModel();
		renderer.setShowEndCarbons(true);
		renderer.setShowExplicitHydrogens(true);
		renderer.setShowImplicitHydrogens(false);
		renderer.setZoomFactor(0.9);

		// Automatically scroll to the center
		JViewport vp = getScrollPane().getViewport();
		vp.setViewPosition(new Point(150, 300));

	}

	/**
	 * Override method to restrict functionality
	 */
	@Override
	public void setupPopupMenus(PopupController2D inputAdapter) {
	}

	/**
	 * Override method to restrict functionality
	 */
	@Override
	public void stateChanged(ChangeEvent e) {

		super.stateChanged(e);

		if (statusLabel == null)
			return;

		// Get the formula and mass of complete molecule
		IChemModel model = getJChemPaintModel().getChemModel();
		IAtomContainer wholeModel = model.getBuilder().newAtomContainer();
		Iterator containers = ChemModelManipulator.getAllAtomContainers(model)
				.iterator();

		while (containers.hasNext()) {
			wholeModel.add((IAtomContainer) containers.next());
		}

		MFAnalyser formulaAnalyzer = new MFAnalyser(wholeModel, true);

		String wholeFormula = formulaAnalyzer
				.getHTMLMolecularFormulaWithCharge();
		double wholeMass = formulaAnalyzer.getMass();

		StringBuilder status = new StringBuilder("<html>Formula: ");
		status.append(wholeFormula);
		status.append(",  mass: ");
		status.append(massFormater.format(wholeMass));
		status.append(" amu");

		Renderer2DModel rendererModel = getJChemPaintModel().getRendererModel();

		IAtomContainer selectedPart = rendererModel.getSelectedPart();

		if ((selectedPart != null) && (selectedPart.getAtomCount() > 0)) {

			MFAnalyser selectionAnalyzer = new MFAnalyser(selectedPart, true);

			String selectionFormula = selectionAnalyzer
					.getHTMLMolecularFormulaWithCharge();
			double selectionMass = selectionAnalyzer.getMass();

			status.append("; selected formula: ");
			status.append(selectionFormula);
			status.append(", mass: ");
			status.append(massFormater.format(selectionMass));
			status.append(" amu");

		}

		status.append("</html>");

		statusLabel.setText(status.toString());

	}

}
