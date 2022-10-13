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
public class SpectrumsViewer {
	public static void main(String[] args) throws SecurityException, IOException, InterruptedException {

		final InputParser setup = new InputParser("spectrums", args);

		// set up the logging
		Logger logger = setUpLogger(setup.filename);
		logger.log(Level.INFO, "beginning "+setup.numRuns +" simulations on "+setup.numCores+" cores");
		logger.log(Level.INFO, "results will be saved to '"+setup.filename+".csv'.");

		final double[] eBins = CSV.readColumn(new File("input/energy.txt"));
		final double[] tBins = CSV.readColumn(new File("input/"+setup.implosionName+" time.txt"));
		final double[][] spectrum = SpectrumGenerator.modifySpectrum(
				SpectrumGenerator.interpretSpectrumFile(
						tBins, eBins,
						CSV.read(new File("input/"+setup.implosionName+" spectrum.txt"), '\t')),
				setup.yieldFactor);
		
		double backgroundExcess = 4e17/Math2.sum(spectrum);

		ExecutorService threads = Executors.newFixedThreadPool(setup.numCores);

		double[][] timeVectors = new double[setup.numRuns][];
		double[][] yieldVectors = new double[setup.numRuns][];
		double[][] temperatureVectors = new double[setup.numRuns][];
		double[][] densityVectors = new double[setup.numRuns][];

		for (int k = 0; k < setup.numRuns; k ++) {
			final int K = k;

			Callable<Void> task = () -> {
				Analysis mc;
				try {
					mc = new Analysis(
							setup.opticsConfig,
							setup.detectorConfig,
							setup.ion,
							setup.shielding*backgroundExcess,
							setup.uncertainty*1e-2,
							false,
							setup.energyBin, setup.timeBin,
							setup.tolerance, logger); // make the simulation
				} catch (IOException e) {
					e.printStackTrace();
					return null;
				}

				logger.log(Level.INFO, String.format("%d/%d",
				                                     K, setup.numRuns));

				try {
					mc.respondAndAnalyze(
							eBins,
							tBins,
							spectrum,
							ErrorMode.HESSIAN); // and run it many times!
				} catch (Exception e) {
					logger.log(Level.SEVERE, e.getMessage(), e);
				}

				timeVectors[K] = mc.getTimeAxis();
				densityVectors[K] = mc.getArealDensity();
				temperatureVectors[K] = mc.getIonTemperature();
				yieldVectors[K] = mc.getNeutronYield();
				if (timeVectors[K] == null) {
					timeVectors[K] = new double[0];
					densityVectors[K] = new double[0];
					temperatureVectors[K] = new double[0];
					yieldVectors[K] = new double[0];
				}

				if (K%10 == 9) {
					try {
						save(timeVectors, yieldVectors, temperatureVectors, densityVectors, setup.filename, logger);
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
		save(timeVectors, yieldVectors, temperatureVectors, densityVectors, setup.filename, logger); // and then, finally, save the result
	}

	private static void save(double[][] timeVectors, double[][] yieldVectors,
	                         double[][] temperatureVectors, double[][] densityVectors,
	                         String filepath, Logger logger) {
		try {
			CSV.write(timeVectors, new File(filepath+"_time.csv"), ',');
			CSV.write(yieldVectors, new File(filepath+"_yield.csv"), ',');
			CSV.write(temperatureVectors, new File(filepath+"_temperature.csv"), ',');
			CSV.write(densityVectors, new File(filepath+"_density.csv"), ',');
			logger.log(Level.INFO, "Saved ensemble results to '"+filepath+"'.");
		} catch (IOException e) {
			logger.log(Level.SEVERE, e.getMessage(), e);
		}
	}
}
