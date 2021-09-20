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
package app;

import java.io.File;
import java.io.IOException;
import java.util.logging.ConsoleHandler;
import java.util.logging.FileHandler;
import java.util.logging.Handler;
import java.util.logging.Level;
import java.util.logging.Logger;

import physics.Analysis;
import physics.Particle;
import physics.SpectrumGenerator;
import util.COSYMapping;
import util.CSV;
import physics.Analysis.ErrorMode;


/**
 * the class that handles the GUI.
 * 
 * @author Justin Kunimune
 */
public class ConsoleEvaluator {
	public static void main(String[] args) throws SecurityException, IOException {
		for (int i = 0; i < Analysis.HEADERS.length; i ++) {
			if (i < 4)
				Analysis.HEADERS_WITH_ERRORS[i] = Analysis.HEADERS[i];
			else {
				Analysis.HEADERS_WITH_ERRORS[2*(i-4)+4] = Analysis.HEADERS[i];
				Analysis.HEADERS_WITH_ERRORS[2*(i-4)+5] = Analysis.HEADERS[i] + " error";
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
		Logger logger = Logger.getLogger("app");
		logger.setUseParentHandlers(false);
		logger.setLevel(Level.ALL);
		Handler consoleHandler = new ConsoleHandler();
		consoleHandler.setLevel(Level.FINE);
		logger.addHandler(consoleHandler);
		Handler logfileHandler = new FileHandler("output/"+filename+".log");
		logger.addHandler(logfileHandler);
		logger.log(Level.INFO, "beginning "+numYields+" evaluations on "+numThreads+" cores");
		
		Thread[] threads = new Thread[numThreads];
		double[][] results = new double[numYields][Analysis.HEADERS_WITH_ERRORS.length];
		for (int t = 0; t < numThreads; t ++) {
			final int T = t;
			final String finalImplosionName = implosionName;
			threads[t] = new Thread(() -> {
				Analysis mc;
				try {
					COSYMapping map = CSV.readCosyCoefficients(new File("input/MRSt_IRF_FP tilted_final.txt"), 3);
					map.setConfig(Particle.D, 12.45);
					mc = new Analysis(
							3e-3,
							2*foilRadius,
							2*foilRadius,
							foilThickness,
							6e0,
							apertureWidth,
							apertureHeight,
							map,
							68,
							.1,
							logger); // make the simulation
				} catch (Exception e) {
					logger.log(Level.SEVERE, e.getMessage(), e);
					return;
				}
				
				for (int k = 0; k < numYields/numThreads; k ++) {
					double[] eBins, tBins;
					double[][] spec;
					try {
						eBins = CSV.readColumn(new File("input/energy.txt"));
						tBins = CSV.readColumn(new File("input/time "+finalImplosionName+".txt"));
						spec = CSV.read(new File("input/spectrum "+finalImplosionName+".txt"), '\t');
						if (spec.length != eBins.length-1 || spec[0].length != tBins.length-1) {
							System.out.println("interpreting a weird spectrum file...");
							spec = SpectrumGenerator.interpretSpectrumFile(tBins, eBins, spec);
						}
					} catch (ArrayIndexOutOfBoundsException | NumberFormatException | IOException e) {
						logger.log(Level.SEVERE, e.getMessage(), e);
						return;
					}

					double yield = Math.pow(10, -3.*Math.random());
					SpectrumGenerator.modifySpectrum(tBins, eBins, spec, yield, 1, 1, 0);
					
					logger.log(Level.INFO, String.format("Yn = %f (%d/%d)", yield,
							T+numThreads*k, numYields));
					
					double[] result;
					try {
						result = mc.respondAndAnalyze(
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
					
					if (T == 0 && (k+1)%5 == 0)
						save(results, filename, logger);
				}
			});
			threads[t].start();
		}
		
		for (int t = 0; t < numThreads; t ++) // wait for all threads to finish
			while (threads[t].isAlive()) {}
		save(results, filename, logger); // and then, finally, save the result
	}
	
	private static void save(double[][] results, String filename, Logger logger) {
		try {
			CSV.write(results, new File("output/"+filename+".csv"), ',', Analysis.HEADERS_WITH_ERRORS);
			logger.log(Level.INFO, "Saved ensemble results to output/"+filename+".csv");
		} catch (IOException e) {
			logger.log(Level.SEVERE, e.getMessage(), e);
		}
	}
}
