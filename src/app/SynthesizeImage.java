/*
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

import physics.Detector.DetectorConfiguration;
import physics.IonOptics;
import physics.IonOptics.IonOpticConfiguration;
import physics.Particle;
import physics.SpectrumGenerator;
import util.Math2;
import util.PythonPlot;

import java.io.IOException;
import java.util.Random;

public class SynthesizeImage {

	public static void main(String[] args) throws IOException {
		IonOpticConfiguration config = IonOpticConfiguration.LOW_RESOLUTION;
		DetectorConfiguration detector = DetectorConfiguration.DOUBLE_STREAK_CAMERA;
		IonOptics optics = new IonOptics(config, Particle.D, detector.cosyFile,
		                                 detector.tiltAngle, detector.offset,
		                                 0.1, false);
		int numSlits = detector.slitPositions.length;

		double[] slitTimes = new double[numSlits];
		for (int s = 0; s < numSlits; s ++) {
			double slitEnergy = optics.energyVsPosition.evaluate(detector.slitPositions[s]);
			slitEnergy /= 1e6*-Particle.E.charge * 8/9.;
			slitTimes[s] = optics.map(slitEnergy)[2];
		}

		Random random = new Random();

		double pixelEdge = 100e-6;//25e-6;
		double resolution = 102e-6;
		double gain = 122.4;
		double[][] xBins = new double[numSlits][];
		for (int s = 0; s < numSlits; s ++) {
			xBins[s] = new double[(int) (detector.slitLengths[0]/pixelEdge) + 1];
			for (int i = 0; i < xBins[s].length; i++) {
				xBins[s][i] = detector.slitLengths[s]*i/(xBins[s].length - 1.) - detector.slitLengths[s]/2 + detector.slitPositions[s];
				xBins[s][i] /= 1e-2;
			}
		}
		double[] yBins = new double[xBins[0].length];
		for (int i = 0; i < yBins.length; i ++) {
			yBins[i] = detector.slitLengths[0]*i/(xBins[0].length - 1.) - detector.slitLengths[0]/2;
			yBins[i] /= 1e-2;
		}
		double[][][] counts = new double[numSlits][][];
		for (int s = 0; s < numSlits; s ++) {
			counts[s] = new double[yBins.length - 1][xBins[s].length - 1];
			for (int i = 0; i < counts[s].length; i ++)
				for (int j = 0; j < counts[s][i].length; j ++)
					counts[s][i][j] = Math2.poisson(100*Math.pow(pixelEdge/25e-6, 2), random);
		}

		double[] energyBins = new double[101];
		for (int i = 0; i < energyBins.length; i ++)
			energyBins[i] = 12 + 4.*i/(energyBins.length - 1);
		double[] spectrum = SpectrumGenerator.generateSpectrum(
				1, 7., 7., 0., .8, energyBins);

		long total = 0, detected = 0;
		for (int k = 0; k < 4e17*optics.efficiency(14); k ++) {
			if (k%10000 == 0)
				System.out.println(k+"/"+(int)(4e17*optics.efficiency(14)));
			double time = Math2.normal(0, 100e-12, random);
			double energy = Math2.drawFromProbabilityDistribution(
					energyBins, spectrum, random);
			double[] position = optics.simulate(energy, time, true);
//			System.out.printf("[%.5f,%.5f],\n", position[0], position[1]);
			for (int s = 0; s < numSlits; s ++) {
				double y = Math2.normal(
					  position[1],
					  resolution, random);
				if (Math.abs(y) < detector.slitWidths[s]/2) {
					double x = Math2.normal(
						  position[0],
						  resolution, random);
					if (Math.abs(x - detector.slitPositions[s]) < detector.slitLengths[s]/2) {
						y += (position[2] - slitTimes[s])/detector.streakTime*detector.slitLengths[s];
						int i = Math2.bin(y/1e-2, yBins);
						int j = Math2.bin(x/1e-2, xBins[s]);
						if (i >= 0 && j >= 0)
							counts[s][i][j] += gain;
						detected ++;
					}
				}
			}
			total ++;
		}
		System.out.printf("out of 4e17 total particles, %d made it thru the ion optics and %d of those were detected (for efficiencies of %.3g and %.3g)\n", total, detected, (float) total/4e17, (float) detected/total);

		for (int s = 0; s < counts.length; s ++)
			PythonPlot.plotHeatmap(xBins[s], yBins, counts[s],
			                       "x (cm)",
			                       "y (cm)",
			                       "Camera "+s+" image");
	}

}
