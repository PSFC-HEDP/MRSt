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
import util.Math2;

import java.io.File;
import java.io.IOError;
import java.io.IOException;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.logging.Level;
import java.util.logging.Logger;

import static util.InputParser.setUpLogger;


/**
 * evaluate error bars at a large number of yields
 * 
 * @author Justin Kunimune
 */
public class ConfigurationEvaluator {
	public static void main(String[] args) throws SecurityException, IOException, InterruptedException {

		final InputParser setup = new InputParser("ensemble", args);

		// set up the logging
		Logger logger = setUpLogger(setup.filename);
		logger.log(Level.INFO, "beginning "+setup.numRuns +" evaluations on "+setup.numCores+" cores");
		logger.log(Level.INFO, "results will be saved to '"+setup.filename+".csv'.");

		final double[] eBins = CSV.readColumn(new File("input/energy.txt"));
		final double[] tBins = CSV.readColumn(new File("input/"+setup.implosionName+" time.txt"));
		final double[][] spectrum = SpectrumGenerator.interpretSpectrumFile(
				tBins, eBins,
				CSV.read(new File("input/"+setup.implosionName+" spectrum.txt"), '\t')
		);

		String[] parameterNames = Analysis.KeyParameterSet.EXAMPLE.getHeaderSansUnits();
		String[] header = new String[1 + 2*parameterNames.length];
		header[0] = "true total yield";
		System.arraycopy(Analysis.KeyParameterSet.EXAMPLE.getHeaderWithUnitsAndError(),
		                 0, header, 1, header.length - 1);
		double[][] results = new double[setup.numRuns][header.length];
		ExecutorService threads = Executors.newFixedThreadPool(setup.numCores);

		for (int k = 0; k < setup.numRuns; k ++) {
			final int K = k;

			Callable<Void> task = () -> {
				double newYield = 1e+19*Math.pow(10, -3.0*Math.random());
				double[][] scaledSpectrum = SpectrumGenerator.modifySpectrum(
						spectrum, newYield/Math2.sum(spectrum));

				Analysis mc;
				try {
					mc = new Analysis(setup, scaledSpectrum, logger); // make the simulation
				} catch (IOException e) {
					e.printStackTrace();
					return null;
				}

				logger.log(Level.INFO, String.format("Yn = %.4g (%d/%d)", newYield,
													 K, setup.numRuns));

				Analysis.KeyParameterSet result;
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

				results[K][0] = newYield;
				if (result != null) {
					for (int i = 0; i < parameterNames.length; i ++) {
						results[K][2*i + 1] = result.getValue(parameterNames[i]);
						results[K][2*i + 2] = Math.sqrt(result.getVariance(parameterNames[i]));
					}
				}
				for (int i = 1; i < results[K].length; i++)
					results[K][i] = Double.NaN;

				if (K % 10 == 9) {
					try {
						save(header, results, setup.filename + ".csv", logger);
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
		save(header, results, setup.filename + ".csv", logger); // and then, finally, save the result
	}
	
	private static void save(String[] header, double[][] results, String filepath, Logger logger) {
		try {
			CSV.write(results, new File(filepath), ',', header);
			logger.log(Level.INFO, "Saved ensemble results to '"+filepath+"'.");
		} catch (IOException e) {
			logger.log(Level.SEVERE, e.getMessage(), e);
		}
	}
}
