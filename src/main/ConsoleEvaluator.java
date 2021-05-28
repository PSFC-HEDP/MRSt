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
	
	
	public static final void main(String[] args) throws SecurityException, IOException {
		for (int i = 0; i < MRSt.HEADERS.length; i ++) {
			if (i < 4)
				MRSt.HEADERS_WITH_ERRORS[i] = MRSt.HEADERS[i];
			else {
				MRSt.HEADERS_WITH_ERRORS[2*(i-4)+4] = MRSt.HEADERS[i];
				int locationOfTheWordQuoteErrorUnquote = MRSt.HEADERS[i].indexOf('(') - 1;
				if (locationOfTheWordQuoteErrorUnquote == -2)
					locationOfTheWordQuoteErrorUnquote = MRSt.HEADERS[i].length();
				MRSt.HEADERS_WITH_ERRORS[2*(i-4)+5] = MRSt.HEADERS[i] + " error";
			}
		}

		double foilRadius = Double.parseDouble(args[0])*1e-6; // foil radius measured in μm
		double foilThickness = Double.parseDouble(args[1])*1e-6; // foil thickness measured in μm
		double apertureWidth = Double.parseDouble(args[2])*1e-3; // aperture width in mm
		double apertureHeight = Double.parseDouble(args[3])*1e-3; // aperture height in mm
		int numYields = Integer.parseInt(args[4]); // number of datums to do
		String implosionName = "og";
		if (args.length > 5)
			implosionName = args[5];
		int numThreads = Math.min(10, Runtime.getRuntime().availableProcessors());
		
		String filename = String.format("ensemble_%.0f_%.0f_%.0f_%.0f_%d_%tF", foilRadius/1e-4, foilThickness/1e-5, apertureWidth/1e-3, apertureHeight/1e-2, numYields, System.currentTimeMillis());
		
		System.setProperty("java.util.logging.SimpleFormatter.format",
				"%1$tF %1$tT | %4$-7s | %5$s%6$s%n");
		Logger logger = Logger.getLogger("main");
		logger.setUseParentHandlers(false);
		logger.setLevel(Level.ALL);
		Handler consoleHandler = new ConsoleHandler();
		consoleHandler.setLevel(Level.FINE);
		logger.addHandler(consoleHandler);
		Handler logfileHandler = new FileHandler("working/"+filename+".log");
		logger.addHandler(logfileHandler);
		logger.log(Level.INFO, "beginning "+numYields+" evaluations on "+numThreads+" cores");
		
		Thread[] threads = new Thread[numThreads];
		double[][] results = new double[numYields][MRSt.HEADERS_WITH_ERRORS.length];
		for (int t = 0; t < numThreads; t ++) {
			final int T = t;
			final String finalImplosionName = implosionName;
			threads[t] = new Thread(() -> {
				MRSt mc = null;
				try {
					COSYMapping map = CSV.readCosyCoefficients(new File("data/MRSt_IRF_FP tilted_final.txt"), 3);
					double[][] cosyCoefficients = map.coefficients;
					int[][] cosyExponents = map.exponents;
					mc = new MRSt(
							Particle.D,
							3e-3,
							2*foilRadius,
							2*foilRadius,
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
				
				for (int k = 0; k < numYields/numThreads; k ++) {
					double[] eBins = null, tBins = null;
					double[][] spec = null;
					try {
						eBins = CSV.readColumn(new File("data/energy.txt"));
						tBins = CSV.readColumn(new File("data/time "+finalImplosionName+".txt"));
						spec = CSV.read(new File("data/spectrum "+finalImplosionName+".txt"), '\t');
						if (spec.length != eBins.length-1 || spec[0].length != tBins.length-1) {
							System.out.println("interpreting a weird spectrum file...");
							spec = MRSt.interpretSpectrumFile(tBins, eBins, spec);
						}
					} catch (ArrayIndexOutOfBoundsException | NumberFormatException e) {
						logger.log(Level.SEVERE, e.getMessage(), e);
					} catch (IOException e) {
						logger.log(Level.SEVERE, e.getMessage(), e);
					}
					
					double yield = Math.pow(10, -3.*Math.random());
					MRSt.modifySpectrum(tBins, eBins, spec, yield, 1, 1, 0);
					
					logger.log(Level.INFO, String.format("Yn = %f (%d/%d)", yield,
							T+numThreads*k, numYields));
					
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
					results[T+numThreads*k][0] = yield;
					results[T+numThreads*k][1] = 1;
					results[T+numThreads*k][2] = 1;
					results[T+numThreads*k][3] = 0;
					if (result != null)
						System.arraycopy(result, 0, results[T+numThreads*k], 4, result.length);
					else
						for (int i = 4; i < results[T+numThreads*k].length; i ++)
							results[T+numThreads*k][i] = Double.NaN;
					
					if ((k+1)%10 == 0) {
//						while (!isAvailable(new File("working/"+filename+".csv"))) {}
						save(results, filename, logger);
					}
				}
			});
			threads[t].start();
		}
		
		for (int t = 0; t < numThreads; t ++)
			while (threads[t].isAlive()) {}
		save(results, filename, logger);
	}
	
	private static void save(double[][] results, String filename, Logger logger) {
		try {
			CSV.write(results, new File("working/"+filename+".csv"), ',', MRSt.HEADERS_WITH_ERRORS);
			logger.log(Level.INFO, "Saved ensemble results to working/"+filename+".csv");
		} catch (IOException e) {
			logger.log(Level.SEVERE, e.getMessage(), e);
		}
	}
}
