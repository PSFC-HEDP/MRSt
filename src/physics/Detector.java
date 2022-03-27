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
public abstract class Detector {

	/**
	 * the background level in a single space unit (signal/space)
	 */
	private final double backgroundPerPixel;
	/**
	 * the background variance in a single space unit (signal^2/space)
	 */
	private final double noisePerPixel;
	/**
	 * the greatest amount of signal in a space unit that is measured linearly (signal/space)
	 */
	private final double saturationLimitPerPixel;

	/**
	 * supply the basic parameters for a generic Detector
	 * @param background the background level in a single space unit (signal/space)
	 * @param noise the background variance in a single space unit (signal^2/space)
	 * @param saturation the signal density at which it saturates (signal/space)
	 */
	protected Detector(double background, double noise, double saturation) {
		this.backgroundPerPixel = background;
		this.noisePerPixel = noise;
		this.saturationLimitPerPixel = saturation;
	}

	/**
	 * compute the detected spectrum given a deuteron spectrum at the photocathode
	 * @param stochastic whether to apply noise to the result
	 */
	abstract double[][] response(double[] energyBins, double[] timeBins,
								 double[][] inSpectrum, boolean stochastic);

	/**
	 * the total of signal units for every incident deuteron
	 * @param energy energy of the neutron (MeV)
	 * @return the gain (signal/deuteron)
	 */
	abstract double gain(double energy);

	/**
	 * the size of one bin in whatever unit is relevant to the detector; this
	 * will be used to scale the background and saturation limit
	 * @param energy the energy at which to evaluate it (MeV)
	 * @param energyBins the energy bin edges (MeV)
	 * @param timeBins the time bin edges (ns)
	 * @return a conversion factor (space/bin)
	 */
	protected double pixelsPerBin(double energy, double[] energyBins, double[] timeBins) {
		return 1;
	}

	/**
	 * the background variance in a single bin given the bins
	 * @param energy the energy at which to evaluate it (MeV)
	 * @param energyBins the energy bin edges (MeV)
	 * @param timeBins the time bin edges (ns)
	 * @return a variance (signal^2/bin)
	 */
	public final double noise(double energy, double[] energyBins, double[] timeBins) {
		return this.noisePerPixel
			  *this.pixelsPerBin(energy, energyBins, timeBins);
	}

	/**
	 * the base signal level in a single bin given the bins
	 * @param energy the energy at which to evaluate it (MeV)
	 * @param energyBins the energy bin edges (MeV)
	 * @param timeBins the time bin edges (ns)
	 */
	public final double background(double energy, double[] energyBins, double[] timeBins) {
		return this.backgroundPerPixel
			  *this.pixelsPerBin(energy, energyBins, timeBins);
	}

	/**
	 * the signal level in one bin that will cause it to start saturating.
	 * @param energyBins the energy bin edges (MeV)
	 * @param timeBins the time bin edges (ns)
	 */
	public final double saturationLimit(double energy, double[] energyBins, double[] timeBins) {
		return this.saturationLimitPerPixel
			  *this.pixelsPerBin(energy, energyBins, timeBins);
	}

	/**
	 * the mean gain over the given energies
	 * @param minEnergy the lower neutron energy bound (MeV n)
	 * @param maxEnergy the upper neutron energy bound (MeV n)
	 */
	public double averageGain(double minEnergy, double maxEnergy) {
		double totalGain = 0;
		int res = 20;
		for (int i = -res; i <= res; i ++) // perform an unweited mean over [13.5, 14.5] MeV
			totalGain += this.gain(minEnergy + (maxEnergy - minEnergy)*i/res);
		return totalGain / (2.*res + 1);
	}

	/**
	 * the fraction of deuterons in this energy range that make it thru
	 * @param minEnergy the lower neutron energy bound (MeV n)
	 * @param maxEnergy the upper neutron energy bound (MeV n)
	 * @return a number in [0, 1]
	 */
	public final double quantumEfficiency(double minEnergy, double maxEnergy) {
		double survivingDeuterons = 0;
		int res = 20;
		for (int i = -res; i <= res; i ++) // perform an unweited mean over [13.5, 14.5] MeV
			if (this.gain(minEnergy + (maxEnergy - minEnergy)*i/res) > 0)
				survivingDeuterons += 1;
		return survivingDeuterons / (2.*res + 1);
	}


	public static class DetectorConfiguration {
		public static DetectorConfiguration SINGLE_STREAK_CAMERA =
			  new DetectorConfiguration("MRSt_IRF_FP_00deg",
										0.00000, 11.5e-9,
										new double[] {-0.75e-2},
										new double[] {2.5e-2},
										new double[] {500e-6});
		public static DetectorConfiguration DOUBLE_STREAK_CAMERA =
			  new DetectorConfiguration("MRSt_IRF_FP_70deg",
										66.58565, 4.5e-9,
										new double[] {-5e-2, 0},
										new double[] {2.5e-2, 2.5e-2},
										new double[] {500e-6, 500e-6});
		public static DetectorConfiguration DOWNSCATTER_SLIT =
			  new DetectorConfiguration("MRSt_IRF_FP_60deg",
										62.44417, 4.5e-9,
										new double[] {-11.5e-2, 0},
										new double[] {2.5e-2, 2.5e-2},
										new double[] {500e-6, 500e-6});
		public static DetectorConfiguration DRIFT_TUBE =
			  new DetectorConfiguration("MRSt_IRF_FP_70deg",
										66.58565, Double.NaN,
										new double[] {0},
										new double[] {1},
										new double[] {1});

		public final String cosyFile;
		public final double tiltAngle;
		public final double streakTime;
		public final double[] slitPositions;
		public final double[] slitLengths;
		public final double[] slitWidths;

		public DetectorConfiguration(
			  String cosyFile, double tiltAngle, double streakTime,
			  double[] slitPositions, double[] slitLengths, double[] slitWidths) {
			this.cosyFile = cosyFile;
			this.tiltAngle = tiltAngle;
			this.streakTime = streakTime;
			this.slitPositions = slitPositions;
			this.slitLengths = slitLengths;
			this.slitWidths = slitWidths;
		}
	}

}
