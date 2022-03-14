package physics;

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
public interface Detector {

	/**
	 * the number of signal electrons created for every incident deuteron
	 * @param energy energy of the neutron [MeV]
	 */
	double efficiency(double energy);

	/**
	 * the average number of cable electrons from each signal electron
	 */
	double gain();

	/**
	 * the signal variance in a single bin given the bins
	 * @param energy the energy at which to evaluate it (MeV)
	 * @param energyBins the energy bin edges (MeV)
	 * @param timeBins the time bin edges (ns)
	 */
	double noise(double energy, double[] energyBins, double[] timeBins);

	/**
	 * the base signal level in a single bin given the bins
	 * @param energy the energy at which to evaluate it (MeV)
	 * @param energyBins the energy bin edges (MeV)
	 * @param timeBins the time bin edges (ns)
	 */
	double background(double energy, double[] energyBins, double[] timeBins);

	/**
	 * compute the detected spectrum given a deuteron spectrum at the photocathode
	 * @param stochastic whether to apply noise to the result
	 */
	double[][] response(double[] energyBins, double[] timeBins,
						double[][] inSpectrum, boolean stochastic);

	class DetectorConfiguration {
		public static DetectorConfiguration SINGLE_STREAK_CAMERA =
			  new DetectorConfiguration( 0.00000, 11.5e-9,
										 new double[] {-0.75e-2},
										 new double[] {2.5e-2},
										 new double[] {500e-6});
		public static DetectorConfiguration DOUBLE_STREAK_CAMERA =
			  new DetectorConfiguration(66.58565,  4.5e-9,
										new double[] {-5e-2, 0},
										new double[] {2.5e-2, 2.5e-2},
										new double[] {500e-6, 500e-6});
		public static DetectorConfiguration DOWNSCATTER_SLIT =
			  new DetectorConfiguration(66.58565,  4.5e-9,
										new double[] {-10e-2, 0},
										new double[] {2.5e-2, 2.5e-2},
										new double[] {500e-6, 500e-6});

		public final double tiltAngle;
		public final double streakTime;
		public final double[] slitPositions;
		public final double[] slitLengths;
		public final double[] slitWidths;

		public DetectorConfiguration(
			  double tiltAngle, double streakTime,
			  double[] slitPositions, double[] slitLengths, double[] slitWidths) {
			this.tiltAngle = tiltAngle;
			this.streakTime = streakTime;
			this.slitPositions = slitPositions;
			this.slitLengths = slitLengths;
			this.slitWidths = slitWidths;
		}
	}

}
