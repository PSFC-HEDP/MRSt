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
package app;

import physics.Analysis;
import physics.Detector.DetectorConfiguration;
import physics.IonOptics.IonOpticConfiguration;

import java.io.IOException;

/**
 * @author Justin Kunimune
 */
public class ConfigurationPlotter {
	public static void main(String[] args) throws NumberFormatException, IOException {
			for (double rFoil = 0; rFoil < 451e-6; rFoil += 500e-6) {
				for (double tFoil = 0; tFoil < 91e-6; tFoil += 15e-6) {
					for (double wAperture = 0; wAperture < 5.1e-3; wAperture += 1e-3) {
//					System.out.println("setting up simulation");

//					double[] errs = new double[MEASUREMENTS.length];
//					double[] eBins = CSV.readColumn(new File("input/energy.txt"));
//					double[] tBins = CSV.readColumn(new File("input/time og.txt"));
//					double[][] spec = CSV.read(new File("input/time og.txt"), '\t');
//					if (spec.length != eBins.length-1 || spec[0].length != tBins.length-1) {
//						System.out.println("interpreting a weird spectrum file...");
//						spec = SpectrumGenerator.interpretSpectrumFile(tBins, eBins, spec);
//					}

					Analysis mc = new Analysis(
							new IonOpticConfiguration(tFoil, rFoil, wAperture),
							DetectorConfiguration.DRIFT_TUBE,
							0,
							false,
							null); // make the simulation

					double[] resolutions = mc.computeResolution(14);
					double timeRes = resolutions[1]; // [ps]
					double energyRes = resolutions[0]; // [keV]

//					for (int i = 0; i < NUM_RUNS; i ++) {
//						double[] result = null;
//						while (result == null) {
//							try {
//								result =mc.respondAndAnalyze(
//															eBins,
//															tBins,
//															spec,
//															ErrorMode.NONE); // and run it many times!
//							}
//							catch (IllegalStateException e) {
//								System.out.println("fuck");
//							}
//						}
//						for (int j = 0; j < errs.length; j ++)
//							errs[j] += Math.pow(result[MEASUREMENTS[j]] - TARGETS[j], 2)/NUM_RUNS;
//					}
//					for (int j = 0; j < errs.length; j ++) {
//						errs[j] = Math.sqrt(errs[j]);
//					}
					
					System.out.printf("[%g, %g, %g, %g, %g, %g, %g],\n",
							rFoil, tFoil, wAperture, 20e-3, energyRes, timeRes, mc.efficiency(14));
				}
			}
		}
	}
}
