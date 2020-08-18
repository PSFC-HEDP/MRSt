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
	
	
	public static final void main(String[] args) {
		char config = args[0].charAt(0);
		int numYields = Integer.parseInt(args[1]);
		if (config != 'h' && config != 'm' && config != 'l')
			throw new IllegalArgumentException("first argument must be 'low', 'med', or 'high'.");
		
		Logger logger = Logger.getLogger("main");
		logger.setLevel(Level.INFO);
//		logger.setLevel(Level.ALL);
		Handler consoleHandler = new ConsoleHandler();
		consoleHandler.setLevel(Level.ALL);
		logger.addHandler(consoleHandler);
		logger.log(Level.INFO, "beginning "+numYields+" evaluations for configuration "+config);
		
		String filename = String.format("ensemble_%s_%d_%tF_%tR", config, numYields, System.currentTimeMillis(), System.currentTimeMillis());
		
		MRSt mc = null;
		try {
			COSYMapping map = CSV.readCosyCoefficients(new File("data/MRSt_IRF_FP tilted_final.txt"), 3);
			double[][] cosyCoefficients = map.coefficients;
			int[][] cosyExponents = map.exponents;
			mc = new MRSt(
					Particle.D,
					3e-3,
					(config == 'h') ? .1e-3 : (config == 'm') ? .2e-3 : .3e-3,
					(config == 'h') ? 25e-6 : (config == 'm') ? 50e-6 : 80e-6,
					CSV.read(new File("data/stopping_power_deuterons.csv"), ','),
					6e0,
					(config == 'h') ? 2e-3 : (config == 'm') ? 3e-3 : 4e-3,
					20e-3,
					COSY_MINIMUM_ENERGY,
					COSY_MAXIMUM_ENERGY,
					COSY_REFERENCE_ENERGY,
					cosyCoefficients,
					cosyExponents,
					66.586,
					(config == 'h') ? 400 : (config == 'm') ? 36 : 1,
					logger); // make the simulation
		} catch (Exception e) {
			logger.log(Level.SEVERE, e.getMessage(), e);
		}
		
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
		
		double[][] results = new double[numYields][HEADERS_WITH_ERRORS.length];
		for (int k = 0; k < numYields; k ++) {
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
			
			if ((k+1)%20 == 0 || k+1 == numYields) {
				try {
					CSV.write(results, new File("working/"+filename), ',', HEADERS_WITH_ERRORS);
					logger.log(Level.INFO, "Saved ensemble results to working/"+filename);
				} catch (IOException e) {
					logger.log(Level.SEVERE, e.getMessage(), e);
				}
			}
		}
	}
}
