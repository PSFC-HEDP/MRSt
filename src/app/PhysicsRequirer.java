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
import java.util.logging.ConsoleHandler;
import java.util.logging.FileHandler;
import java.util.logging.Formatter;
import java.util.logging.Handler;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.logging.SimpleFormatter;


/**
 * attempt to discern a small number of implosions from one another
 *
 * @author Justin Kunimune
 */
public class PhysicsRequirer {
	public static void main(String[] args) throws SecurityException, IOException, InterruptedException {

		// first, parse the program arguments
		final InputParser setup = new InputParser("comparison", args);

		if (setup.opticsConfig == null)
			throw new IllegalArgumentException("you need to always specify the ion optic configuration from now on.");
		if (setup.detectorConfig == null)
			throw new IllegalArgumentException("you need to always specify the detector configuration from now on.");

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
		logger.log(Level.INFO, "beginning "+setup.numYields+" evaluations on "+setup.numCores+" cores");
		logger.log(Level.INFO, "results will be saved to '"+setup.filename+".csv'.");

		final double[] eBins = CSV.readColumn(new File("input/scan/energy.txt"));
		final double[] tBins = CSV.readColumn(new File("input/scan/time.txt"));

		ExecutorService threads = Executors.newFixedThreadPool(setup.numCores);

		String[] headers = new String[Analysis.HEADERS_WITH_ERRORS.length + 1];
		headers[0] = "Case";
		System.arraycopy(Analysis.HEADERS_WITH_ERRORS, 0, headers, 0, Analysis.HEADERS_WITH_ERRORS.length);
		List<double[]> results = new ArrayList<double[]>();

		int caseIndex = 0;
		int runIndex = 0;
		for (Path path : (Iterable<Path>)Files.walk(Paths.get("input/scan/"))::iterator) {
			if (path.getFileName().toString().startsWith("spectrum ")) {
				final int J = caseIndex;

				String key = path.getFileName().toString();
				key = key.substring(9, key.length() - 4);

				double[][] initialSpec = CSV.read(new File("input/spectrum " + key + ".txt"), '\t');
				final double[][] spectrum;
				if (initialSpec.length != eBins.length - 1 || initialSpec[0].length != tBins.length - 1) {
					System.out.println("interpreting a weird spectrum file...");
					spectrum = SpectrumGenerator.interpretSpectrumFile(tBins, eBins, initialSpec);
				}
				else {
					spectrum = initialSpec;
				}

				for (int dk = 0; dk < setup.numYields; dk++) {
					final int K = runIndex;
					final boolean saveHere = dk == setup.numYields - 1;

					threads.submit(() -> {
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
							return;
						}

						double yield = 1e+19*Math.pow(10, -3.0*Math.random());
						double[][] scaledSpectrum = SpectrumGenerator.modifySpectrum(
							  spectrum, yield);

						logger.log(Level.INFO, String.format("Yn = %.4g (%d/%d)", yield,
															 K, setup.numYields));

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
							double[] resultSlot = results.get(K);
							resultSlot[0] = J;
							resultSlot[1] = yield;
							if (result != null)
								System.arraycopy(result, 0, resultSlot, 2, result.length);
							else
								for (int i = 1; i < results.get(K).length; i++)
									results.get(K)[i] = Double.NaN;
						} catch (IndexOutOfBoundsException e) {
							logger.log(Level.SEVERE, e.getMessage(), e);
						}
						if (saveHere) {
							try {
								save(results, setup.filename + ".csv", headers, logger);
							} catch (IOError e) {
								logger.log(Level.SEVERE, e.getMessage(), e);
							}
						}
					});
				}
				caseIndex ++;
			}
		}
		threads.shutdown();
		threads.awaitTermination(3, TimeUnit.DAYS); // wait for all threads to finish
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
