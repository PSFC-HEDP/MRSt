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
import physics.Detector.DetectorConfiguration;
import physics.IonOptics.IonOpticConfiguration;
import physics.Particle;
import physics.SpectrumGenerator;
import util.COSYMapping;
import util.CSV;

import java.io.File;
import java.io.IOException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.logging.ConsoleHandler;
import java.util.logging.FileHandler;
import java.util.logging.Formatter;
import java.util.logging.Handler;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.logging.SimpleFormatter;


/**
 * the class that handles the GUI.
 * 
 * @author Justin Kunimune
 */
public class ConsoleEvaluator {
	public static void main(String[] args) throws SecurityException, IOException, InterruptedException {
		StringBuilder filename = new StringBuilder("ensemble");
		String prospectiveImplosionName = "og with falling temp";
		int prospectiveNumYields = 1000;
		int numCores = Math.min(10, Runtime.getRuntime().availableProcessors());
		IonOpticConfiguration prospectiveIonConfig = IonOpticConfiguration.HIGH_EFFICIENCY;
		DetectorConfiguration prospectiveDetectorConfig = DetectorConfiguration.SINGLE_STREAK_CAMERA;
		double prospectiveUncertainty = 0;
		for (String arg : args) {
			if (arg.contains("=")) {
				String key = arg.substring(0, arg.indexOf('='));
				String value = arg.substring(arg.indexOf('=') + 1);
				switch (key) {
					case "implosion":
						prospectiveImplosionName = value;
						break;
					case "yields":
						prospectiveNumYields = Integer.parseInt(value);
						break;
					case "cores":
						numCores = Integer.parseInt(value);
						break;
					case "uncertainty":
						prospectiveUncertainty = Double.parseDouble(value);
						break;
					case "optics":
						if (value.toLowerCase().startsWith("h"))
							prospectiveIonConfig = IonOpticConfiguration.HIGH_EFFICIENCY;
						else if (value.toLowerCase().startsWith("m"))
							prospectiveIonConfig = IonOpticConfiguration.MID_EFFICIENCY;
						else if (value.toLowerCase().startsWith("l"))
							prospectiveIonConfig = IonOpticConfiguration.LOW_EFFICIENCY;
						else
							System.err.println("I don't know the '" + value + "' configuration");
						break;
					case "detector":
						if (value.toLowerCase().startsWith("1"))
							prospectiveDetectorConfig = DetectorConfiguration.SINGLE_STREAK_CAMERA;
						else if (value.toLowerCase().startsWith("2"))
							prospectiveDetectorConfig = DetectorConfiguration.DOUBLE_STREAK_CAMERA;
						else if (value.toLowerCase().startsWith("d"))
							prospectiveDetectorConfig = DetectorConfiguration.DOWNSCATTER_SLIT;
						else
							System.err.println("I don't know the '" + value + "' camera");
						break;
					default:
						System.err.println("I don't know to what '" + key + "' refers");
						break;
				}

				String tag = value.replace(" ", "");
				if (tag.length() > 6) tag = tag.substring(0, 6);
				filename.append("_").append(tag);
			}
			else {
				System.err.println("I don't understand '"+arg+"'.");
			}
		}

		final String implosionName = prospectiveImplosionName;
		final int numYields = prospectiveNumYields;
		final DetectorConfiguration detectorConfiguration = prospectiveDetectorConfig;
		final IonOpticConfiguration ionOpticConfiguration = prospectiveIonConfig;
		final double uncertainty = prospectiveUncertainty;
		
		filename.append(String.format("_%tF", System.currentTimeMillis()));

		// set up the logging
		System.setProperty("java.util.logging.SimpleFormatter.format",
						   "%1$tF %1$tT | %4$-7s | %5$s%6$s%n");
		// (I don't know why they make this so difficult)
		Formatter formatter = new SimpleFormatter();
		Logger logger = Logger.getLogger("app");
		logger.setUseParentHandlers(false);
		logger.setLevel(Level.ALL);
		Handler consoleHandler = new ConsoleHandler();
		consoleHandler.setLevel(Level.FINE);
		consoleHandler.setFormatter(formatter);
		logger.addHandler(consoleHandler);
		Handler logfileHandler = new FileHandler("output/"+filename+".log");
		logfileHandler.setLevel(Level.INFO);
		logfileHandler.setFormatter(formatter);
		logger.addHandler(logfileHandler);
		logger.log(Level.INFO, "beginning "+prospectiveNumYields+" evaluations on "+numCores+" cores");

		final COSYMapping map;
		if (detectorConfiguration.tiltAngle == 0)
			map = CSV.readCosyCoefficients(new File("input/MRSt_IRF_FP not tilted.txt"), 3);
		else
			map = CSV.readCosyCoefficients(new File("input/MRSt_IRF_FP tilted.txt"), 3);
		map.setConfig(Particle.D, 12.45);

		final double[] eBins = CSV.readColumn(new File("input/energy.txt"));
		final double[] tBins = CSV.readColumn(new File("input/time "+implosionName+".txt"));
		double[][] initialSpec = CSV.read(new File("input/spectrum "+implosionName+".txt"), '\t');
		final double[][] spectrum;
		if (initialSpec.length != eBins.length-1 || initialSpec[0].length != tBins.length-1) {
			System.out.println("interpreting a weird spectrum file...");
			spectrum = SpectrumGenerator.interpretSpectrumFile(tBins, eBins, initialSpec);
		}
		else
			spectrum = initialSpec;

		ExecutorService threads = Executors.newFixedThreadPool(numCores);
		double[][] results = new double[prospectiveNumYields][Analysis.HEADERS_WITH_ERRORS.length];
		for (int k = 0; k < numYields; k ++) {
			final int K = k;

			threads.submit(() -> {
				Analysis mc;
				try {
					mc = new Analysis(
						  ionOpticConfiguration,
						  detectorConfiguration,
						  uncertainty*1e-2,
						  false,
						  .1,
						  logger); // make the simulation
				} catch (IOException e) {
					e.printStackTrace();
					return;
				}

				double yield = 1e+19*Math.pow(10, -3.0*Math.random());
				double[][] scaledSpectrum = SpectrumGenerator.modifySpectrum(
					  spectrum, yield);

				logger.log(Level.INFO, String.format("Yn = %.4g (%d/%d)", yield,
						K, numYields));

				double[] result;
				try {
					result = mc.respondAndAnalyze(
							eBins,
							tBins,
							scaledSpectrum,
							ErrorMode.HESSIAN); // and run it many times!
				} catch (Exception e) {
					logger.log(Level.SEVERE, e.getMessage(), e);
					result = null;
				}
				results[K][0] = yield;
				if (result != null)
					System.arraycopy(result, 0, results[K], 1, result.length);
				else
					for (int i = 1; i < results[K].length; i ++)
						results[K][i] = Double.NaN;

				if (K%10 == 9)
					save(results, filename.toString(), logger);
			});
		}

		threads.shutdown();
		threads.awaitTermination(3, TimeUnit.DAYS); // wait for all threads to finish
		save(results, filename.toString(), logger); // and then, finally, save the result
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
