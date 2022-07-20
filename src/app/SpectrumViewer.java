/**
 * MIT License
 * 
 * Copyright (c) 2018 Justin Kunimune
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
package app;

import physics.Analysis;
import physics.Analysis.ErrorMode;
import physics.Detector.DetectorConfiguration;
import physics.IonOptics.IonOpticConfiguration;
import physics.SpectrumGenerator;
import util.CSV;
import util.Math2;
import util.PythonPlot;

import java.io.File;
import java.io.IOException;
import java.util.logging.ConsoleHandler;
import java.util.logging.Level;
import java.util.logging.Logger;


/**
 * the class that handles the GUI.
 * 
 * @author Justin Kunimune
 */
public class SpectrumViewer {
	/**
	 * build the GUI and display it.
	 */
	public static void main(String[] args) throws NumberFormatException, IOException {
		IonOpticConfiguration optics = IonOpticConfiguration.HIGH_EFFICIENCY;
		DetectorConfiguration detector = DetectorConfiguration.DRIFT_TUBE;
		String trajectory = "scan/p2 p1 p4 burnoff";
		double yieldFactor = 1;
		boolean reuseMatrix = false;

		System.setProperty("java.util.logging.SimpleFormatter.format",
		                   "%1$tF %1$tT | %4$-7s | %5$s%6$s%n");
		Logger logger = Logger.getLogger("app");
		logger.setUseParentHandlers(false);
		logger.setLevel(Level.ALL);
		ConsoleHandler commandlineHandler = new ConsoleHandler();
		commandlineHandler.setLevel(Level.ALL);
		logger.addHandler(commandlineHandler);

		double[] energyBins = CSV.readColumn(new File("input/energy.txt"));
		double[] timeBins = CSV.readColumn(new File("input/" + trajectory + " time.txt"));
		double[][] spectrum = CSV.read(new File("input/" + trajectory + " spectrum.txt"), '\t');
		ErrorMode errorMode = ErrorMode.HESSIAN;

		double[] eBins, tBins;
		double[][] spec;
		try {
			eBins = energyBins.clone(); // save the current values of these spectra
			tBins = timeBins.clone();
			spec = deepClone(spectrum);
			spec = SpectrumGenerator.interpretSpectrumFile(tBins, eBins, spec);
			spec = SpectrumGenerator.modifySpectrum(spec, yieldFactor*Math2.sum(spec));
		} catch (ArrayIndexOutOfBoundsException e) {
			logger.severe("Invalid input spectrum file.");
			return;
		}

		logger.log(Level.INFO, "running fit on spectrum with yield factor = "+yieldFactor);
		Analysis mc;
		try {
			mc = new Analysis(
					optics,
					detector,
					0,
					reuseMatrix,
					logger); // make the simulation

			double dispersion = mc.computeDispersion();
			double skew = mc.computeTimeSkew();
			logger.info(String.format("Dispersion: %.2f keV/mm", dispersion));
			logger.info(String.format("Time skew:  %.2f ps/keV", skew));
//			double[] res = mc.computeResolution(14.);
//			logger.info(String.format("Energy res: %.2f keV", res[0]));
//			logger.info(String.format("Time res:   %.2f ps", res[1]));

			mc.respondAndAnalyze(
					eBins,
					tBins,
					spec,
					errorMode); // and run it!

		} catch (Exception e) {
			logger.log(Level.SEVERE, e.getMessage(), e);
			return;
		}

		double[][] smallSpec = Math2.downsample(tBins, eBins, spec, mc.getTimeBins(), mc.getEnergyBins());
		try { // send the data to python for plotting
			PythonPlot.plotHeatmap(mc.getTimeBins(), mc.getEnergyBins(), smallSpec,
			                       "Time (ns)", "Energy (MeV)", "Original neutron spectrum");
			PythonPlot.plotHeatmap(mc.getTimeBins(), mc.getDeuteronEnergyBins(), mc.getSignalDistribution(),
			            "Time (ns)", "Energy (MeV)", "Synthetic signal distribution");
//			PythonPlot.plotHeatmap(mc.getTimeBins(), mc.getEnergyBins(), mc.getFitNeutronSpectrum(),
//			            "Fitted neutron spectrum");
//			PythonPlot.plotHeatmap(mc.getTimeBins(), mc.getEnergyBins(), mc.getFitSignalDistribution(),
//			            "Fitted signal distribution");
			PythonPlot.plotLines("Trajectories", trajectory + "{}",
					mc.getTimeAxis(), "Time (ns)",
			        mc.getNeutronYield(), mc.getNeutronYieldError(), "Yn (10^15/ns)",
					mc.getIonTemperature(), mc.getIonTemperatureError(), "Ti (keV)"//,
//					mc.getArealDensity(), mc.getArealDensityError(), "ρR (g/cm^2)"
			);
			PythonPlot.compareHeatmap(mc.getTimeBins(), mc.getDeuteronEnergyBins(), mc.getSignalDistribution(), mc.getFitSignalDistribution(),
					"Time (ns)", "Energy (MeV)", "Synthetic signal", "Fit signal");
//			PythonPlot.compareHeatmap(mc.getTimeBins(), mc.getDeuteronEnergyBins(), mc.getDeuteronSpectrum(), mc.getFitDeuteronSpectrum(),
//									  "Time", "Energy (MeV)", "Synthetic deuteron spectrum", "Fit deuteron spectrum");
			PythonPlot.compareHeatmap(mc.getTimeBins(), mc.getEnergyBins(), smallSpec, mc.getFitNeutronSpectrum(),
									  "Time (ns)", "Energy (MeV)", "Original neutron spectrum", "Fit neutron spectrum");
			PythonPlot.compareHeatmap(mc.getTimeBins(), mc.getDeuteronEnergyBins(),
			                          mc.efficiencyCorrect(mc.getSignalDistribution(), true),
			                          mc.efficiencyCorrect(mc.getFitSignalDistribution(), false),
			                          mc.efficiencyCorrect(mc.getIdealSignalDistribution(), false),
			                          "Time (ns)", "Energy (MeV)",
			                          "Corrected signal distribution", "Fit signal distribution", "True spectrum");
			PythonPlot.compareHeatmap(mc.getTimeBins(), mc.getDeuteronEnergyBins(),
									  mc.getSignalDistribution(), mc.getBackgroundSpectrum(),
									  "Time (ns)", "Energy (MeV)", "Synthetic signal distribution", "Detector background");
		} catch (IOException e) {
			logger.log(Level.SEVERE, "Could not access plotting scripts and/or plots", e);
		}
	}

	private static double[][] deepClone(double[][] old) {
		double[][] now = new double[old.length][old[0].length];
		for (int i = 0; i < old.length; i ++)
			for (int j = 0; j < old[i].length; j ++)
				now[i][j] = old[i][j];
		return now;
	}
	
}
