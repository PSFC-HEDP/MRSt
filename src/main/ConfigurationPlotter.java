/**
 * MIT License
 * 
 * Copyright (c) 2018 Justin Kunimune
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
import main.CSV.COSYMapping;
import main.Analysis.ErrorMode;

/**
 * @author Justin Kunimune
 */
public class ConfigurationPlotter {
	
	private static final int[] MEASUREMENTS = {6, 22, 24}; // indices for burn width, vi, and dTdt
	private static final double[] TARGETS = {.066711, 67.33, 1.83};
	private static final int NUM_RUNS = 4;
	private static final int NUM_PARTICLES = 10000;
	private static final int RESOLUTION = 100;
	private static final double COSY_MINIMUM_ENERGY = 10.7e6;
	private static final double COSY_MAXIMUM_ENERGY = 14.2e6;
	private static final double COSY_REFERENCE_ENERGY = 12.45e6;
	
	
	/**
	 * @param args
	 * @throws IOException 
	 * @throws NumberFormatException 
	 */
	public static void main(String[] args) throws NumberFormatException, IOException {
		for (double rFoil = 200e-6; rFoil < 201e-6; rFoil += 100e-6) {
			for (double tFoil = 25e-6; tFoil < 110e-6; tFoil += 15e-6) {
				for (double wAperture = 1.0e-3; wAperture < 5.1e-3; wAperture += 1.0e-3) {
//					System.out.println("setting up simulatin");
					Analysis mc = null;
					COSYMapping map = CSV.readCosyCoefficients(new File("data/MRSt_IRF_FP tilted_final.txt"), 3);
					double[][] cosyCoefficients = map.coefficients;
					int[][] cosyExponents = map.exponents;
					mc = new Analysis(
							Particle.D,
							3e-3,
							2*rFoil,
							2*rFoil,
							tFoil,
							CSV.read(new File("data/stopping_power_deuterons.csv"), ','),
							6e0,
							wAperture,
							20e-3,
							COSY_MINIMUM_ENERGY,
							COSY_MAXIMUM_ENERGY,
							COSY_REFERENCE_ENERGY,
							cosyCoefficients,
							cosyExponents,
							68,
							.1,
							null); // make the simulation
					
//					System.out.println("kalkula emablia");

					double[] resolutions = mc.computeResolution(14);
					double timeRes = resolutions[1]; // [ps]
					double energyRes = resolutions[0]; // [keV]
					
//					System.out.println("kalkula veria");
					double[] errs = new double[MEASUREMENTS.length];
					double[] eBins = null, tBins = null;
					double[][] spec = null;
					eBins = CSV.readColumn(new File("data/Energy bins.txt"));
					tBins = CSV.readColumn(new File("data/nsp_150327_16p26_time - copia.txt"));
					spec = CSV.read(new File("data/nsp_150327_16p26.txt"), '\t');
					if (spec.length != eBins.length-1 || spec[0].length != tBins.length-1) {
						System.out.println("interpreting a weird spectrum file...");
						spec = SpectrumGenerator.interpretSpectrumFile(tBins, eBins, spec);
					}
					
					for (int i = 0; i < NUM_RUNS; i ++) {
						double[] result = null;
						while (result == null) {
							try {
								result =mc.respondAndAnalyze(
															eBins,
															tBins,
															spec,
															ErrorMode.NONE); // and run it many times!
							}
							catch (IllegalStateException e) {
								System.out.println("fuck");
							}
						}
						for (int j = 0; j < errs.length; j ++)
							errs[j] += Math.pow(result[MEASUREMENTS[j]] - TARGETS[j], 2)/NUM_RUNS;
					}
					for (int j = 0; j < errs.length; j ++) {
						errs[j] = Math.sqrt(errs[j]);
					}
					
					System.out.printf("[%g, %g, %g, %g, %g, %g, %g, %g, %g, %g],\n",
							rFoil, tFoil, wAperture, 20e-3, energyRes, timeRes, mc.efficiency(), errs[0], errs[1], errs[2]);
				}
			}
		}
	}
}
