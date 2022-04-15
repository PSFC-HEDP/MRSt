/**
 * MIT License
 * <p>
 * Copyright (c) 2021 Justin Kunimune
 * <p>
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * <p>
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * <p>
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
import physics.Particle;
import physics.SpectrumGenerator;
import util.COSYMapping;
import util.CSV;
import util.Math2;
import util.PythonPlot;

import java.io.File;
import java.io.IOException;

public class FocalPlanePlotter {
	public static void main(String[] args) throws IOException {
		// select constants
		DetectorConfiguration slits = DetectorConfiguration.DRIFT_TUBE;
		COSYMapping cosyMapping = CSV.readCosyCoefficients(
			  new File(String.format("input/%s.txt", slits.cosyFile)),
			  3, Particle.D, 12.45);

		// set up the simulation
		IonOptics io = new IonOptics(
				3e-3,
				.8e-3,
				.8e-3,
				90e-6,
				6e0,
				5e-3,
				20e-3,
				12,
				16,
				cosyMapping,
				slits.tiltAngle,
				0,
				false
		);

		// select the energies
		int N = 35;
		double[] energiesN = new double[N + 1];
		double[] energiesD = new double[N + 1];
		for (int i = 0; i <= N; i ++) {
			energiesD[i] = 10.7 + 3.5*i/N;
			energiesN[i] = energiesD[i]*9/8.;
		}
		energiesN[0] = 12;
		energiesN[N] = 16;

		// then collect the values
		int M = 1000;
		double[][] positions = new double[N + 1][3*M];
		for (int i = 0; i <= N; i ++) {
			System.out.printf("%d/%d: %.2f MeV <- %.2f MeV\n", i, N + 1, energiesD[i], energiesN[i]);
			for (int k = 0; k < M; k ++) {
				double[] xyzt = io.simulate(energiesN[i], 0, false);
				positions[i][3*k  ] = xyzt[0]/Math.cos(Math.toRadians(slits.tiltAngle));
				positions[i][3*k+1] = xyzt[1];
				positions[i][3*k+2] = xyzt[3];
			}
		}

		// now come up with a simple spectrum to illustrate the particle densities
		double[] kernel = new double[15];
		for (int i = 0; i < kernel.length; i ++) {
			int di = i - kernel.length/2;
			double dEdi = energiesN[1] - energiesN[0];
			kernel[i] = Math.exp(-Math.pow((di*dEdi + .200)/.200, 2)/2.);
		}
		double[] spectrum = Math2.convolve(
				SpectrumGenerator.generateSpectrum(
						1, 12., 12., 0, 1.0, energiesN
				),
				kernel
		);

		// finally, plot
		PythonPlot.plotFocalPlane(
			  energiesD, positions, spectrum,
			  slits.slitPositions, slits.slitLengths, slits.slitWidths);
	}
}
