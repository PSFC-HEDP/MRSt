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

import util.Math2;

import static physics.Analysis.NOISE_RANDOM;

public class RealisticDetector extends Detector {
	private final double timeResolution;


	protected RealisticDetector(IonOptics optics, double shielding, double timeResolution) {
		// background number comes from Wink 2016
		super(optics.backgroundAtDetector()/shielding,
		      optics.backgroundAtDetector()/shielding,
		      Double.POSITIVE_INFINITY, 1);
		this.timeResolution = timeResolution;
	}

	@Override
	public double[][] response(double[] energyBins, double[] timeBins, double[][] inSpectrum,
	                           boolean stochastic, boolean background, boolean gaps) {
		// compute time resolution degradation
		double timeStep = timeBins[1] - timeBins[0];
		double[] timeResponse;
		if (timeBins.length%2 == 1)
			timeResponse = new double[timeBins.length];
		else
			timeResponse = new double[timeBins.length - 1];
		for (int j = 0; j < timeResponse.length; j ++) {
			int dj = j - timeResponse.length/2;
			if (timeResolution >= timeStep/10.)
				timeResponse[j] = Math.exp(-Math.pow(dj*timeStep/(timeResolution/2.355), 2)/2.);
			else
				timeResponse[j] = (dj == 0) ? 1 : 0;
		}
		double timeTotal = Math2.sum(timeResponse);
		for (int j = 0; j < timeResponse.length; j ++)
			timeResponse[j] /= timeTotal;

		// apply time resolution degradation
		double[][] sansBackground = new double[inSpectrum.length][inSpectrum[0].length];
		for (int i = 0; i < inSpectrum.length; i ++) {
			for (int l = 0; l < timeResponse.length; l ++) {
				int dj = l - timeResponse.length/2;
				for (int j = 0; j < timeBins.length - 1; j ++) {
					if (j + dj >= 0 && j + dj < timeBins.length - 1) {
						if (stochastic)
							sansBackground[i][j] += gain*Math2.binomial(
									(int)inSpectrum[i][j + dj], timeResponse[l], NOISE_RANDOM); // as well as gain
						else
							sansBackground[i][j] += gain*inSpectrum[i][j + dj]*timeResponse[l];
					}
				}
			}
		}
		// add background
		double[][] withBackground = new double[inSpectrum.length][inSpectrum[0].length];
		for (int i = 0; i < inSpectrum.length; i ++) {
			double energy = (energyBins[i] + energyBins[i + 1])/2;
			double level = background(energy, energyBins, timeBins, gaps);
			for (int j = 0; j < inSpectrum[i].length; j ++) {
				if (stochastic)
					withBackground[i][j] = sansBackground[i][j] + Math2.poisson(level, NOISE_RANDOM);
				else
					withBackground[i][j] = sansBackground[i][j] + level;
			}
		}
		return withBackground;
	}

	@Override
	double efficiency(double energy, boolean gaps) {
		return 1;
	}

}
