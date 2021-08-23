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

import physics.Analysis;
import physics.IonOptics;
import physics.Particle;
import util.CSV;
import util.CSV.COSYMapping;

import java.io.File;
import java.io.IOException;

public class FocalPlanePlotter {
	public static void main(String[] args) throws IOException {
		COSYMapping map = CSV.readCosyCoefficients(new File("input/MRSt_IRF_FP tilted_final.txt"), 3);
		double[][] cosyCoefficients = map.coefficients;
		int[][] cosyExponents = map.exponents;
		IonOptics io = new IonOptics(
				Particle.D,
				3e-3,
				.8e-3,
				.8e-3,
				90e-6,
				6e0,
				5e-3,
				20e-3,
				12,
				16,
				14,
				cosyCoefficients,
				cosyExponents,
				68
		); // make the simulation
		System.out.print("Ys = np.array(\n");
		for (double E = 12; E <= 16; E += .125) {
			System.out.println("],[");
			for (int k = 0; k < 1000; k++) {
				io.simulate(E, 0);
			}
		}
		System.out.print("])\n");
		System.out.print("E = np.array([");
		for (double E = 12; E <= 16; E += .125) {
			System.out.printf("%.3f, ", E);
		}
		System.out.print("])\n");
		io.simulate(13.54, 0);
		io.simulate(14.46, 0);
	}
}
