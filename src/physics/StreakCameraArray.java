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
package physics;

import util.COSYMapping;
import util.NumericalMethods.DiscreteFunction;

import java.util.Arrays;

public class StreakCameraArray implements Detector {

	private static final int FP_RESOLUTION = 30;

	private final double slitWidth; // (m)
	private final double[] slitLeftBounds; // (MeV neutron)
	private final double[] slitRiteBounds; // (MeV neutron)

	private final DiscreteFunction bowtieHite; // (MeV neutron -> m)

	public StreakCameraArray(
		  double slitLength, double slitWidth, double slitSpacing,
		  IonOptics optics) {
		this.slitWidth = slitWidth;

		double[] ERef = new double[FP_RESOLUTION];
		double[] xRef = new double[FP_RESOLUTION];
		double[] hRef = new double[FP_RESOLUTION];
		for (int i = 0; i < FP_RESOLUTION; i ++) {
			ERef[i] = 12. + 4.*i/(FP_RESOLUTION - 1); // (MeV neutron)
			double num = 0, xSum = 0, yMax = 0;
			for (int k = 0; k < 1000; k ++) {
				double[] rt = optics.simulate(ERef[i], 0, false);
				if (!Double.isNaN(rt[0])) {
					num += 1;
					xSum += Math.hypot(rt[0], rt[2])*Math.signum(rt[0]);
					if (rt[1] > yMax) yMax = rt[1];
				}
			}
			xRef[i] = xSum/num;
			hRef[i] = yMax*2;
		}
		DiscreteFunction fpPosition = new DiscreteFunction(ERef, xRef, true);
		DiscreteFunction fpEnergy = fpPosition.inv().indexed(FP_RESOLUTION);
		this.bowtieHite = new DiscreteFunction(ERef, hRef, true);

		double fpCenter = fpPosition.evaluate(14);
		double[] slitPositions = {fpCenter - slitSpacing, fpCenter, fpCenter + slitSpacing};
		this.slitLeftBounds = new double[3];
		this.slitRiteBounds = new double[3];
		for (int i = 0; i < this.slitLeftBounds.length; i ++) {
			this.slitLeftBounds[i] = fpEnergy.evaluate(
				  slitPositions[i] - slitLength/2);
			this.slitRiteBounds[i] = fpEnergy.evaluate(
				  slitPositions[i] + slitLength/2);
		}
	}

	@Override
	public double efficiency(double energy) {
		for (int i = 0; i < slitLeftBounds.length; i ++)
			if (energy >= slitLeftBounds[i] && energy <= slitRiteBounds[i])
				return Math.min(1, slitWidth/bowtieHite.evaluate(energy));
		return 0;
	}

	@Override
	public double gain() {
		return 1;
	}

	@Override
	public double[][] response(double[] energyBins, double[] timeBins,
							   double[][] inSpectrum, boolean stochastic) {
		double[][] outSpectrum = new double[energyBins.length-1][timeBins.length-1];
		for (int i = 0; i < energyBins.length-1; i ++) {
			for (int j = 0; j < timeBins.length-1; j ++) {
//				System.out.println((energyBins[i] + energyBins[i+1])/2);
				outSpectrum[i][j] = inSpectrum[i][j]*this.efficiency((energyBins[i] + energyBins[i+1])/2); // TODO: account for time blur from slit width
			}
		}
		return outSpectrum;
	}
}
