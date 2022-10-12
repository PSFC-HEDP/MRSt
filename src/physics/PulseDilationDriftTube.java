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

public class PulseDilationDriftTube extends Detector {
	protected PulseDilationDriftTube(DetectorConfiguration config) {
		// background number comes from Wink 2016
		super(6e5/3.5e16*4e17/config.shielding,
		      6e5/3.5e16*4e17/config.shielding,
		      Double.POSITIVE_INFINITY, 1);
	}

	@Override
	public double[][] response(double[] energyBins, double[] timeBins, double[][] sansBackground,
	                           boolean stochastic, boolean background, boolean gaps) {
//		if (backgroundYield < 4e16 || backgroundYield > 5e18)
//			System.out.println(backgroundYield);
		double[][] withBackground = new double[sansBackground.length][];
		for (int i = 0; i < sansBackground.length; i ++) {
			withBackground[i] = new double[sansBackground[i].length];
			double energy = (energyBins[i] + energyBins[i + 1])/2;
			double level = background(energy, energyBins, timeBins, gaps);
			for (int j = 0; j < sansBackground[i].length; j ++) {
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
