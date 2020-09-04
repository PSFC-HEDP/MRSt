/**
 * MIT License
 * 
 * Copyright (c) 2020 Justin Kunimune
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
package main;

import java.io.File;
import java.io.IOException;
import java.util.logging.ConsoleHandler;
import java.util.logging.FileHandler;
import java.util.logging.Handler;
import java.util.logging.Level;
import java.util.logging.Logger;
import main.CSV.COSYMapping;
import main.MRSt.ErrorMode;


/**
 * the class that handles the GUI.
 * 
 * @author Justin Kunimune
 */
public class ConsoleEvaluator {
	
	private static final double COSY_MINIMUM_ENERGY = 10.7e6;
	private static final double COSY_MAXIMUM_ENERGY = 14.2e6;
	private static final double COSY_REFERENCE_ENERGY = 12.45e6;
	
	private static final String[] HEADERS = {
		"Yield factor", "Temperature factor", "Down-scatter factor", "Velocity shift (μm/ns)",
		"Computation time (s)", "Total yield (10^15)", "Bang time (ns)",
		"Burn width (ns)", "Burn skew", "Burn kurtosis",
		"Max \u03C1R (ns)", "\u03C1R at max (g/cm^2)",
		"Burn-average \u03C1R (g/cm^2)", "d\u03C1R/dt at BT (g/cm^2/ns)",
		"Burn-average Ti (keV)", "dTi/dt at BT (keV/ns)",
		"Burn-average vi (km/s)", "dvi/dt at BT (km/s/ns)"};
	private static final String[] HEADERS_WITH_ERRORS = new String[(HEADERS.length-4)*2+4];
	
	
	public static final void main(String[] args) throws SecurityException, IOException {
		for (int i = 0; i < HEADERS.length; i ++) {
			if (i < 4)
				HEADERS_WITH_ERRORS[i] = HEADERS[i];
			else {
				HEADERS_WITH_ERRORS[2*(i-4)+4] = HEADERS[i];
				int locationOfTheWordQuoteErrorUnquote = HEADERS[i].indexOf('(') - 1;
				if (locationOfTheWordQuoteErrorUnquote == -2)
					locationOfTheWordQuoteErrorUnquote = HEADERS[i].length();
				HEADERS_WITH_ERRORS[2*(i-4)+5] = HEADERS[i] + " error";
			}
		}

		double foilRadius = Double.parseDouble(args[0])*1e-6; // foil radius measured in μm
		double foilThickness = Double.parseDouble(args[1])*1e-6; // foil thickness measured in μm
		double apertureWidth = Double.parseDouble(args[2])*1e-3; // aperture width in mm
		double apertureHeight = Double.parseDouble(args[3])*1e-3; // aperture height in mm
		int numYields = Integer.parseInt(args[4]); // number of datums to do
		
		String filename = String.format("ensemble_%.0f_%.0f_%.1f_%.0f_%d_%tF", foilRadius/1e-4, foilThickness/1e-5, apertureWidth/1e-3, apertureHeight/1e-2, numYields, System.currentTimeMillis());
		
		Logger logger = Logger.getLogger("main");
		logger.setUseParentHandlers(false);
		logger.setLevel(Level.ALL);
		Handler consoleHandler = new ConsoleHandler();
		consoleHandler.setLevel(Level.FINER);
		logger.addHandler(consoleHandler);
		Handler logfileHandler = new FileHandler("working/"+filename+".log");
		logger.addHandler(logfileHandler);
		logger.log(Level.INFO, "beginning "+numYields+" evaluations");
		
		MRSt mc = null;
		try {
			COSYMapping map = CSV.readCosyCoefficients(new File("data/MRSt_IRF_FP tilted_final.txt"), 3);
			double[][] cosyCoefficients = map.coefficients;
			int[][] cosyExponents = map.exponents;
			mc = new MRSt(
					Particle.D,
					3e-3,
					foilRadius,
					foilThickness,
					CSV.read(new File("data/stopping_power_deuterons.csv"), ','),
					6e0,
					apertureWidth,
					apertureHeight,
					COSY_MINIMUM_ENERGY,
					COSY_MAXIMUM_ENERGY,
					COSY_REFERENCE_ENERGY,
					cosyCoefficients,
					cosyExponents,
					68,
					.1,
					logger); // make the simulation
		} catch (Exception e) {
			logger.log(Level.SEVERE, e.getMessage(), e);
		}
		
		double[][] results = new double[numYields][HEADERS_WITH_ERRORS.length];
		for (int k = 0; k < numYields; k ++) {
			double[] eBins = null, tBins = null;
			double[][] spec = null;
			try {
				eBins = CSV.readColumn(new File("data/Energy bins.txt"));
				tBins = CSV.readColumn(new File("data/nsp_150327_16p26_time - copia.txt"));
				spec = CSV.read(new File("data/nsp_150327_16p26.txt"), '\t');
				spec = MRSt.interpretSpectrumFile(tBins, eBins, spec);
			} catch (ArrayIndexOutOfBoundsException | NumberFormatException e) {
				logger.log(Level.SEVERE, e.getMessage(), e);
			} catch (IOException e) {
				logger.log(Level.SEVERE, e.getMessage(), e);
			}
			
			double yield = Math.pow(10, -3.*Math.random());
			MRSt.modifySpectrum(tBins, eBins, spec, yield, 1, 1, 0);
			
			logger.log(Level.INFO, String.format("Yn = %f (%d/%d)", yield, k, numYields));
			
			double[] result;
			try {
				result = mc.respond(
						eBins,
						tBins,
						spec,
						ErrorMode.HESSIAN); // and run it many times!
			} catch (Exception e) {
				logger.log(Level.SEVERE, e.getMessage(), e);
				result = null;
			}
			results[k][0] = yield;
			results[k][1] = 1;
			results[k][2] = 1;
			results[k][3] = 0;
			if (result != null)
				System.arraycopy(result, 0, results[k], 4, result.length);
			else
				for (int i = 4; i < results[k].length; i ++)
					results[k][i] = Double.NaN;
			
			if ((k+1)%5 == 0 || k+1 == numYields) {
				try {
					CSV.write(results, new File("working/"+filename+".csv"), ',', HEADERS_WITH_ERRORS);
					logger.log(Level.INFO, "Saved ensemble results to working/"+filename+".csv");
				} catch (IOException e) {
					logger.log(Level.SEVERE, e.getMessage(), e);
				}
			}
		}
	}
}
