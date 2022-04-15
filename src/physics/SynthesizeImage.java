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
package physics;

import physics.Detector.DetectorConfiguration;
import physics.IonOptics.IonOpticConfiguration;
import util.Math2;
import util.PythonPlot;

import java.io.IOException;
import java.util.Random;

public class SynthesizeImage {

	public static void main(String[] args) throws IOException { // TODO: background
		IonOpticConfiguration config = IonOpticConfiguration.MID_EFFICIENCY;
		DetectorConfiguration detector = DetectorConfiguration.DOWNSCATTER_SLIT;
		IonOptics optics = new IonOptics(config, detector.cosyFile, detector.tiltAngle, 0.1, false);
		int numSlits = detector.slitPositions.length;

		double[] slitTimes = new double[numSlits];
		for (int s = 0; s < numSlits; s ++) {
			double slitEnergy = optics.energyVsPosition.evaluate(detector.slitPositions[s]);
			slitEnergy /= 1e6*-Particle.E.charge * 8/9.;
			slitTimes[s] = optics.map(slitEnergy)[3];
		}

		double pixelEdge = 100e-6;//25e-6;
		double resolution = 102e-6;
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
//		double[] tBins = new double[yBins.length];
//		for (int i = 0; i < yBins.length; i ++)
//			tBins[i] = yBins[i]/detector.slitLengths[0]*detector.streakTime;
		double[][][] counts = new double[numSlits][][];
		for (int s = 0; s < numSlits; s ++)
			counts[s] = new double[yBins.length - 1][xBins[s].length - 1];

		Random random = new Random();
		for (int k = 0; k < 4e17*optics.efficiency(14); k ++) {
			double time = Math2.normal(0, 100e-12, random);
			double energy;
			if (Math.random() < .95)
				energy = Math2.normal(14, .300, random);
			else
				energy = 14 - Math2.gamma(2, 1.0, random);
			double[] position = optics.simulate(energy, time, true);
			for (int s = 0; s < numSlits; s ++) {
				double y = Math2.normal(
					  position[1],
					  resolution, random);
				if (Math.abs(y) < detector.slitWidths[s]/2) {
					double x = Math2.normal(
						  position[0]/Math.cos(Math.toRadians(detector.tiltAngle)),
						  resolution, random);
					if (Math.abs(x - detector.slitPositions[s]) < detector.slitLengths[s]/2) {
						y += (position[3] - slitTimes[s])/detector.streakTime*detector.slitLengths[s];
						int i = Math2.bin(x/1e-2, xBins[s]);
						int j = Math2.bin(y/1e-2, yBins);
//						System.out.println(i + " " + x + "  " + Arrays.toString(xBins[s]));
//						System.out.println(s+","+j + " from " + y + " in  " + Arrays.toString(yBins));
						if (i >= 0 && j >= 0)
							counts[s][j][i]++;
					}
				}
			}
		}

		for (int s = 0; s < counts.length; s ++)
			PythonPlot.plotHeatmap(xBins[s], yBins, counts[s], "Streak camera slit "+s+" image");
	}

}
