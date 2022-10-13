/**
 * MIT License
 *
 * Copyright (c) 2022 Justin Kunimune
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
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.logging.Level;
import java.util.logging.Logger;

import static util.InputParser.setUpLogger;


/**
 * attempt to discern a small number of implosions from one another
 *
 * @author Justin Kunimune
 */
public class PhysicsRequirer {
	public static void main(String[] args) throws SecurityException, IOException, InterruptedException {

		// first, parse the program arguments
		final InputParser setup = new InputParser("comparison", args);

		// set up the logging
		Logger logger = setUpLogger(setup.filename);
		logger.log(Level.INFO, "beginning "+setup.numRuns +" evaluations on "+setup.numCores+" cores");
		logger.log(Level.INFO, "results will be saved to '"+setup.filename+".csv'.");

		final double[] eBins = CSV.readColumn(new File("input/energy.txt"));

		ExecutorService threads = Executors.newFixedThreadPool(setup.numCores);
		
		String[] parameterNames = Analysis.KeyParameterSet.EXAMPLE.getHeaderSansUnits();
		String[] header = new String[1 + 2*parameterNames.length];
		header[0] = "case";
		System.arraycopy(Analysis.KeyParameterSet.EXAMPLE.getHeaderWithUnitsAndError(),
		                 0, header, 1, header.length - 1);
		List<double[]> results = new ArrayList<double[]>();

		for (int runIndex = 0; runIndex < setup.numRuns; runIndex ++) {
			final int K = runIndex;

			int caseIndex = 0;
			for (Path path : (Iterable<Path>)Files.walk(Paths.get("input/"))::iterator) {
				if (path.getFileName().toString().contains("scan") &&
				    path.getFileName().toString().contains("trajectories")) {
					if (runIndex == 0)
						logger.info("Loading scenario " + caseIndex + ": " + path);
					final int J = caseIndex;

					String key = path.getFileName().toString();
					key = key.substring(0, key.length() - 17);

					double[] tBins = CSV.readColumn(new File("input/" + key + " time.txt"));
					double[][] spectrum = SpectrumGenerator.interpretSpectrumFile(
							tBins, eBins,
							CSV.read(new File("input/" + key + " spectrum.txt"), '\t')
					);
					
					double backgroundExcess = 4e17/Math2.sum(spectrum);
					logger.info(String.format("I'm assuming the background levels are all based on 4×10^17, and increasing " +
					                          "the shielding by a factor of %.3g accordingly.",
					                          backgroundExcess));

					threads.submit(() -> {
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
							logger.log(Level.SEVERE, e.getMessage(), e);
							return;
						}

						logger.log(Level.INFO, String.format("Spectrum %d, run %d/%d",
															 J, K, setup.numRuns));

						Analysis.KeyParameterSet result;
						try {
							result = mc.respondAndAnalyze(
								  eBins,
								  tBins,
								  spectrum,
								  ErrorMode.HESSIAN); // and run it many times!
						} catch (Exception e) {
							logger.log(Level.SEVERE, e.getMessage(), e);
							result = null;
						}

						try {
							double[] resultVector = new double[header.length];
							resultVector[0] = J;
							if (result != null) {
								for (int i = 0; i < parameterNames.length; i ++) {
									resultVector[2*i + 1] = result.getValue(parameterNames[i]);
									resultVector[2*i + 2] = Math.sqrt(result.getVariance(parameterNames[i]));
								}
							}
							else
								for (int i = 1; i < resultVector.length; i++)
									resultVector[i] = Double.NaN;
							results.add(resultVector);
						} catch (IndexOutOfBoundsException e) {
							logger.log(Level.SEVERE, e.getMessage(), e);
						}
						if (K%10 == 9) {
							try {
								save(results, setup.filename + ".csv", header, logger);
							} catch (IOError e) {
								logger.log(Level.SEVERE, e.getMessage(), e);
							}
						}
					});

					caseIndex ++;
				}
			}

			if (caseIndex == 0)
				logger.severe("No scenarios found");
		}

		threads.shutdown();
		threads.awaitTermination(3, TimeUnit.DAYS); // wait for all threads to finish
		save(results, setup.filename + ".csv", header, logger); // then save the final result
	}

	private static void save(List<double[]> results, String filepath, String[] headers, Logger logger) {
		try {
			CSV.write(results.toArray(new double[0][]), new File(filepath),
					  ',', headers);
			logger.log(Level.INFO, "Saved ensemble results to '"+filepath+"'.");
		} catch (IOException e) {
			logger.log(Level.SEVERE, e.getMessage(), e);
		}
	}

}
