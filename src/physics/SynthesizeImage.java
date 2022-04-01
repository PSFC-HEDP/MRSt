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

	public static void main(String[] args) throws IOException {
		IonOpticConfiguration config = IonOpticConfiguration.MID_EFFICIENCY;
		DetectorConfiguration detector = DetectorConfiguration.SINGLE_STREAK_CAMERA;
		IonOptics optics = new IonOptics(config, detector.cosyFile, detector.tiltAngle, 0.1, false);

		double pixelEdge = 25e-6;
		double[] xBins = new double[(int)(detector.slitLengths[0]/pixelEdge) + 1];
		double[] yBins = new double[xBins.length];
		for (int i = 0; i < xBins.length; i ++) {
			xBins[i] = detector.slitLengths[0]*i/(xBins.length - 1.) - detector.slitLengths[0]/2 + detector.slitPositions[0];
			xBins[i] /= 1e-2;
		}
		for (int i = 0; i < yBins.length; i ++) {
			yBins[i] = detector.slitLengths[0]*i/(xBins.length - 1.) - detector.slitLengths[0]/2;
			yBins[i] /= 1e-2;
		}
//		double[] tBins = new double[yBins.length];
//		for (int i = 0; i < yBins.length; i ++)
//			tBins[i] = yBins[i]/detector.slitLengths[0]*detector.streakTime;
		double[][] counts = new double[xBins.length - 1][yBins.length - 1];

		double time0 = optics.simulate(13.8, 0, false)[3];
		Random random = new Random();
		for (int k = 0; k < 4e17*optics.efficiency(14); k ++) {
			double time = Math2.normal(0, 100e-12, random);
			double energy;
			if (Math.random() < .95)
				energy = Math2.normal(14, .300, random);
			else
				energy = 14 - Math2.gamma(2, 3.0, random);
			double[] position = optics.simulate(energy, 0, true);
			if (Math.abs(position[1]) < detector.slitWidths[0]/2) {
				double x = position[0]/Math.cos(detector.tiltAngle);
				if (Math.abs(x - detector.slitPositions[0]) < detector.slitLengths[0]/2) {
					double y = position[1] + (position[3] - time0)/detector.streakTime*detector.slitLengths[0];
					int i = Math2.bin(x/1e-2, xBins), j = Math2.bin(y/1e-2, yBins);
//					System.out.println(x+"  "+Arrays.toString(xBins));
//					System.out.println(y+"  "+Arrays.toString(yBins));
					if (i >= 0 && j >= 0)
						counts[j][i] ++;
				}
			}
		}

		PythonPlot.plotHeatmap(yBins, yBins, counts, "Streak camera image");
	}

}
