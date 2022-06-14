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
import physics.SpectrumGenerator;
import util.CSV;
import util.InputParser;

import java.io.File;
import java.io.IOError;
import java.io.IOException;
import java.util.concurrent.Callable;
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
 * evaluate error bars at a large number of yields
 * 
 * @author Justin Kunimune
 */
public class ConsoleEvaluator {
	public static void main(String[] args) throws SecurityException, IOException, InterruptedException {

		final InputParser setup = new InputParser("ensemble", args);

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
		Handler logfileHandler = new FileHandler(setup.filename+".log");
		logfileHandler.setLevel(Level.INFO);
		logfileHandler.setFormatter(formatter);
		logger.addHandler(logfileHandler);
		logger.log(Level.INFO, "beginning "+setup.numRuns +" evaluations on "+setup.numCores+" cores");
		logger.log(Level.INFO, "results will be saved to '"+setup.filename+".csv'.");

		final double[] eBins = CSV.readColumn(new File("input/energy.txt"));
		final double[] tBins = CSV.readColumn(new File("input/time "+setup.implosionName+".txt"));
		final double[][] spectrum = SpectrumGenerator.interpretSpectrumFile(
				tBins, eBins,
				CSV.read(new File("input/spectrum "+setup.implosionName+".txt"), '\t')
		);

		ExecutorService threads = Executors.newFixedThreadPool(setup.numCores);
		double[][] results = new double[setup.numRuns][Analysis.HEADERS_WITH_ERRORS.length];
		for (int k = 0; k < setup.numRuns; k ++) {
			final int K = k;

			Callable<Void> task = () -> {
				Analysis mc;
				try {
					mc = new Analysis(
						  setup.opticsConfig,
						  setup.detectorConfig,
						  setup.uncertainty*1e-2,
						  false,
						  setup.energyBin, setup.timeBin, setup.tolerance,
						  logger); // make the simulation
				} catch (IOException e) {
					e.printStackTrace();
					return null;
				}

				double yield = 1e+19*Math.pow(10, -3.0*Math.random());
				double[][] scaledSpectrum = SpectrumGenerator.modifySpectrum(
					  spectrum, yield);

				logger.log(Level.INFO, String.format("Yn = %.4g (%d/%d)", yield,
													 K, setup.numRuns));

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

				try {
					results[K][0] = yield;
					if (result != null)
						System.arraycopy(result, 0, results[K], 1, result.length);
					else
						for (int i = 1; i < results[K].length; i++)
							results[K][i] = Double.NaN;
				} catch (IndexOutOfBoundsException e) {
					logger.log(Level.SEVERE, e.getMessage(), e);
				}

				if (K%10 == 9) {
					try {
						save(results, setup.filename + ".csv", logger);
					} catch (IOError e) {
						logger.log(Level.SEVERE, e.getMessage(), e);
					}
				}
				return null;
			};
			threads.submit(task);
		}

		threads.shutdown();
		threads.awaitTermination(3, TimeUnit.DAYS); // wait for all threads to finish
		save(results, setup.filename + ".csv", logger); // and then, finally, save the result
	}
	
	private static void save(double[][] results, String filepath, Logger logger) {
		try {
			CSV.write(results, new File(filepath), ',', Analysis.HEADERS_WITH_ERRORS);
			logger.log(Level.INFO, "Saved ensemble results to '"+filepath+"'.");
		} catch (IOException e) {
			logger.log(Level.SEVERE, e.getMessage(), e);
		}
	}
}
