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
import physics.Particle;
import physics.SpectrumGenerator;
import util.CSV;
import util.Math2;
import util.PythonPlot;

import java.io.File;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;

import static util.InputParser.setUpLogger;


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
		Particle ion = Particle.D;
		IonOpticConfiguration optics = IonOpticConfiguration.MEDIUM_RESOLUTION;
		DetectorConfiguration detector = DetectorConfiguration.DRIFT_TUBE;
		String simulationName = "haan shockmerge";
		double yieldFactor = 1;
		boolean reuseMatrix = false;

		Logger logger = setUpLogger(null);

		double[] energyBins = CSV.readColumn(new File("input/energy.txt"));
		double[] timeBins = CSV.readColumn(new File("input/" + simulationName + " time.txt"));
		double[][] spectrum = CSV.read(new File("input/" + simulationName + " spectrum.txt"), '\t');
		ErrorMode errorMode = ErrorMode.HESSIAN;

		double[] eBins, tBins;
		double[][] spec;
		try {
			eBins = energyBins.clone(); // save the current values of these spectra
			tBins = timeBins.clone();
			spec = deepClone(spectrum);
			spec = SpectrumGenerator.interpretSpectrumFile(tBins, eBins, spec);
			spec = SpectrumGenerator.modifySpectrum(spec, yieldFactor);
		} catch (ArrayIndexOutOfBoundsException e) {
			logger.severe("Invalid input spectrum file.");
			return;
		}

		logger.log(Level.INFO, String.format("running fit on spectrum with yield = %.4g × %.4g = %.4g",
		                                     yieldFactor, Math2.sum(spec)/yieldFactor, Math2.sum(spec)));
		Analysis mc;
		try {
			mc = new Analysis(
					optics,
					detector,
					ion,
					293*4e17/Math2.sum(spec),
					0,
					reuseMatrix,
					logger); // make the simulation

			double dispersion = mc.computeDispersion();
			logger.info(String.format("Dispersion: %.2f keV/mm", dispersion));
			double skew = mc.computeTimeSkew();
			logger.info(String.format("Time skew:  %.2f ps/keV", skew));
			double efficiency = mc.efficiency(14);
			logger.info(String.format("Efficiency: %.3g", efficiency));
			double[] res = mc.computeResolution(14.);
			logger.info(String.format("Energy res: %.2f keV", res[0]));
			logger.info(String.format("Time res:   %.2f ps", res[1]));

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
			PythonPlot.plotHeatmap(mc.getDilatedTimeBins(), mc.getDeuteronEnergyBins(), mc.getSignalDistribution(),
			            "Time (ns)", "Energy (MeV)", "Synthetic signal distribution");
			PythonPlot.plotHeatmap(mc.getTimeBins(), mc.getEnergyBins(), mc.getFitNeutronSpectrum(),
			            "Time (ns)", "Energy (MeV)", "Fitted neutron spectrum");
			PythonPlot.plotHeatmap(mc.getDilatedTimeBins(), mc.getDeuteronEnergyBins(), mc.getFitSignalDistribution(),
			            "Time (ns)", "Energy (MeV)", "Fitted signal distribution");
			PythonPlot.plotLines("Trajectories", "input/" + simulationName + " {}.csv",
					mc.getTimeAxis(), "Time (ns)",
			        mc.getNeutronYield(), mc.getNeutronYieldError(), "Yn (10^15/ns)",
					mc.getIonTemperature(), mc.getIonTemperatureError(), "Ti (keV)",
					mc.getArealDensity(), mc.getArealDensityError(), "ρR (g/cm^2)"
			);
			PythonPlot.compareHeatmap(mc.getDilatedTimeBins(), mc.getDeuteronEnergyBins(), mc.getSignalDistribution(), mc.getFitSignalDistribution(),
					"Time (ns)", "Energy (MeV)", "Synthetic signal", "Fit signal");
//			PythonPlot.compareHeatmap(mc.getTimeBins(), mc.getDeuteronEnergyBins(), mc.getDeuteronSpectrum(), mc.getFitDeuteronSpectrum(),
//									  "Time", "Energy (MeV)", "Synthetic deuteron spectrum", "Fit deuteron spectrum");
			PythonPlot.compareHeatmap(mc.getTimeBins(), mc.getEnergyBins(), smallSpec, mc.getFitNeutronSpectrum(),
									  "Time (ns)", "Energy (MeV)", "Original neutron spectrum", "Fit neutron spectrum");
			PythonPlot.compareHeatmap(mc.getDilatedTimeBins(), mc.getDeuteronEnergyBins(),
			                          mc.efficiencyCorrect(mc.getSignalDistribution(), true),
			                          mc.efficiencyCorrect(mc.getFitSignalDistribution(), false),
			                          mc.efficiencyCorrect(mc.getIdealSignalDistribution(), false),
			                          "Time (ns)", "Energy (MeV)",
			                          "Corrected signal distribution", "Fit signal distribution", "True spectrum");
			PythonPlot.compareHeatmap(mc.getDilatedTimeBins(), mc.getDeuteronEnergyBins(),
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
