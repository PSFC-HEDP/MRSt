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

import physics.Analysis;
import physics.Analysis.ErrorMode;
import physics.Particle;
import physics.SpectrumGenerator;
import util.COSYMapping;
import util.CSV;

import java.io.File;
import java.io.IOException;
import java.util.logging.ConsoleHandler;
import java.util.logging.FileHandler;
import java.util.logging.Handler;
import java.util.logging.Level;
import java.util.logging.Logger;


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

		String implosionName;
		if (args.length > 0)
			implosionName = args[0];
		else
			implosionName = "og with falling temp";
		int numYields; // number of datums to do
		if (args.length > 1)
			numYields = Integer.parseInt(args[1]);
		else
			numYields = 1000;
		double tiltAngle;
		if (args.length > 2)
			tiltAngle = Double.parseDouble(args[2]);
		else
			tiltAngle = 66.59;
		double foilRadius, foilThickness, apertureWidth, apertureHeight;
		if (args.length > 6) {
			foilRadius = Double.parseDouble(args[3])*1e-6; // foil radius measured in μm
			foilThickness = Double.parseDouble(args[4])*1e-6; // foil thickness measured in μm
			apertureWidth = Double.parseDouble(args[5])*1e-3; // aperture width in mm
			apertureHeight = Double.parseDouble(args[6])*1e-3; // aperture height in mm
		}
		else {
			if (args.length <= 3 || args[3].equals("h")) {
				foilRadius = 400e-6;
				foilThickness = 90e-6;
				apertureWidth = 5e-3;
				apertureHeight = 20e-3;
			}
			else if (args[3].equals("m")) {
				foilRadius = 200e-6;
				foilThickness = 50e-6;
				apertureWidth = 4e-3;
				apertureHeight = 20e-3;
			}
			else if (args[3].equals("l")) {
				foilRadius = 100e-6;
				foilThickness = 25e-6;
				apertureWidth = 2e-3;
				apertureHeight = 20e-3;
			}
			else {
				throw new IllegalArgumentException(args[3]);
			}
		}
		int numThreads = Math.min(10, Runtime.getRuntime().availableProcessors());
		
		String filename = String.format("ensemble_%.0f_%.0f_%d_%s_%tF",
										apertureWidth/1e-3, tiltAngle, numYields,
										implosionName.replace(" ",""),
										System.currentTimeMillis());
		
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
					COSYMapping map;
					if (tiltAngle == 0)
						map = CSV.readCosyCoefficients(new File("input/MRSt_IRF_FP not tilted.txt"), 3);
					else
						map = CSV.readCosyCoefficients(new File("input/MRSt_IRF_FP tilted.txt"), 3);
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
							tiltAngle,
							false,
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

					double yield = 4e+17*Math.pow(10, -3.0*Math.random());
					SpectrumGenerator.modifySpectrum(spec, yield);
					
					logger.log(Level.INFO, String.format("Yn = %.4g (%d/%d)", yield,
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
					if (result != null)
						System.arraycopy(result, 0, results[T+numThreads*k], 1, result.length);
					else
						for (int i = 1; i < results[T+numThreads*k].length; i ++)
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
