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

/**
 * A detector that perfectly preserves all deuteron informacion, but loses all
 * time resolucion
 */
public class SolidState extends Detector {

	public SolidState() {
		super(0,
			  0,
			  1_000,
			  1);
	}

	@Override
	public double efficiency(double energy) {
		return 1;
	}

	@Override
	public double pixelsPerBin(double energy, double[] energyBins, double[] timeBins) {
		double binSize = (energyBins[1] - energyBins[0]); // (MeV/bin)
		double dispersion = 1; // (m^2/MeV) (TODO)
		double frameArea = 1.36e-7; // m^2/frame (TODO)
		return binSize*dispersion/frameArea;
	}

	@Override
	public double[][] response(double[] energyBins, double[] timeBins,
							   double[][] timeResolved, boolean stochastic,
							   boolean background, boolean gaps) {
		double[][] timeIntegrated = new double[energyBins.length - 1][timeBins.length - 1];
		int j0 = (timeBins.length - 1)/2;
		for (int i = 0; i < timeResolved.length; i ++)
			for (int j = 0; j < timeResolved[i].length; j ++)
				timeIntegrated[i][j0] += timeResolved[i][j];
		return timeIntegrated;
	}
}
