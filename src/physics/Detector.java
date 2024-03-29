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
	 * the background level in a single space unit (signal/space/10^15neutron)
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
	 * the amount of signal generated for each deuteron that is detected (signal/deuteron)
	 */
	public final double gain;

	/**
	 * supply the basic parameters for a generic Detector
	 * @param background the background level in a single space unit (signal/(MeV*ns))
	 * @param noise the background variance in a single space unit (signal^2/(MeV*ns))
	 * @param saturation the signal density at which it saturates (signal/(MeV*ns))
	 * @param gain the amount of signal generated for each deuteron that is detected (signal/deuteron)
	 */
	protected Detector(double background, double noise, double saturation, double gain) {
		this.backgroundPerPixel = background;
		this.noisePerPixel = noise;
		this.saturationLimitPerPixel = saturation;
		this.gain = gain;
	}

	/**
	 * compute the detected spectrum given a deuteron spectrum at the photocathode
	 * @param stochastic whether to apply noise to the result
	 */
	abstract double[][] response(double[] energyBins, double[] timeBins,
								 double[][] inSpectrum, boolean stochastic,
								 boolean background, boolean gaps);

	/**
	 * the total of signal units for every incident deuteron
	 * @param energy energy of the neutron (MeV)
	 * @param gaps whether to set the efficiency to 0 where there are gaps
	 * @return the gain (signal/deuteron)
	 */
	abstract double efficiency(double energy, boolean gaps);

	/**
	 * the size of one bin in whatever unit is relevant to the detector; this
	 * will be used to scale the background and saturation limit
	 * @param energy the energy at which to evaluate it (MeV)
	 * @param energyBins the energy bin edges (MeV)
	 * @param timeBins the time bin edges (ns)
	 * @return a conversion factor (space/bin^2)
	 */
	protected double pixelsPerBin(double energy, double[] energyBins, double[] timeBins) {
		return (energyBins[1] - energyBins[0]) * (timeBins[1] - timeBins[0]);
	}

	/**
	 * the background variance in a single bin given the bins
	 * @param energy the energy at which to evaluate it (MeV)
	 * @param energyBins the energy bin edges (MeV)
	 * @param timeBins the time bin edges (ns)
	 * @return a variance (signal^2/bin^2)
	 */
	public final double noise(double energy, double[] energyBins, double[] timeBins, boolean gaps) {
		if (this.efficiency(energy, gaps) > 0)
			return this.noisePerPixel
				  *this.pixelsPerBin(energy, energyBins, timeBins);
		else
			return 0;
	}

	/**
	 * the base signal level in a single bin given the bins
	 * @param energy the energy at which to evaluate it (MeV)
	 * @param energyBins the energy bin edges (MeV)
	 * @param timeBins the time bin edges (ns)
	 * @return the bin value (signal/bin^2)
	 */
	public final double background(double energy, double[] energyBins, double[] timeBins, boolean gaps) {
		if (this.efficiency(energy, gaps) > 0)
			return this.backgroundPerPixel
				  *this.pixelsPerBin(energy, energyBins, timeBins);
		else
			return 0;
	}

	/**
	 * the signal level in one bin that will cause it to start saturating.
	 * @param energyBins the energy bin edges (MeV)
	 * @param timeBins the time bin edges (ns)
	 * @return the maximum value at which the response is linear (signal/bin^2)
	 */
	public final double saturationLimit(double energy, double[] energyBins, double[] timeBins) {
		return this.saturationLimitPerPixel
			  *this.pixelsPerBin(energy, energyBins, timeBins);
	}


	public static Detector newDetector(DetectorConfiguration config, IonOptics optics, double shielding) {
		if (!Double.isNaN(config.streakTime))
			return new StreakCameraArray(config, optics, shielding);
		else
			return new RealisticDetector(optics, shielding, config.timeResolution);
	}


	public static class DetectorConfiguration {
		public static DetectorConfiguration MAXIMUM_COVERAGE =
			  new DetectorConfiguration("MRSt_IRF_FP_00deg",
										Double.NaN, 0.00000, 0, 11.5e-9,
										new double[] {-7.5e-2, 0.2e-2},
										new double[] {2.5e-2, 2.5e-2},
										new double[] {400e-6, 400e-6});
		public static DetectorConfiguration SINGLE_STREAK_CAMERA =
				new DetectorConfiguration("MRSt_IRF_FP_00deg",
				                          Double.NaN, 0.00000, 0, 11.5e-9,
				                          new double[] {-1.2e-2},
				                          new double[] {2.5e-2},
				                          new double[] {400e-6});
		public static DetectorConfiguration JOHAN =
				new DetectorConfiguration("MRSt_IRF_FP_00deg",
				                          Double.NaN, 0.00000, 0, 5.0e-9,
				                          new double[] {0.0e-2},
				                          new double[] {1.1e-2},
				                          new double[] {400e-6});
		public static DetectorConfiguration DOUBLE_STREAK_CAMERA =
				new DetectorConfiguration("MRSt_IRF_FP_70deg",
				                          Double.NaN, 66.58565, 0, 4.5e-9,
				                          new double[] {-5.5e-2, 1.5e-2},
				                          new double[] {2.5e-2, 2.5e-2},
				                          new double[] {400e-6, 400e-6});
		public static DetectorConfiguration LONG_LONG =
				new DetectorConfiguration("MRSt_IRF_FP_70deg",
				                          Double.NaN, 66.58565, 0, 7.0e-9,
				                          new double[] {-5.5e-2, 1.5e-2},
				                          new double[] {4.0e-2, 4.0e-2},
				                          new double[] {400e-6, 400e-6});
		public static DetectorConfiguration DRIFT_TUBE =
				new DetectorConfiguration("MRSt_IRF_FP_70deg",
				                          66.72555, 66.58565, 0, Double.NaN,
				                          new double[] {},
				                          new double[] {},
				                          new double[] {});
		public static DetectorConfiguration OMEGA =
				new DetectorConfiguration("MRSt_IRF_OMEGA",
				                          63.36917, 63.36917, 0, Double.NaN,
				                          new double[] {},
				                          new double[] {},
				                          new double[] {});
		public static DetectorConfiguration CLEMENT =
				new DetectorConfiguration("MRSt_IRF_FP_70deg",
				                          66.72555, 66.58565, 50e-3, Double.NaN,
				                          new double[] {},
				                          new double[] {},
				                          new double[] {});
		public static DetectorConfiguration PERFECT =
				new DetectorConfiguration("perfect lens",
				                          Double.NaN, 0, 0, Double.NaN,
				                          new double[] {},
				                          new double[] {},
				                          new double[] {});

		public final String cosyFile;
		public final double tiltAngleProton;
		public final double tiltAngleDeuteron;
		public final double streakTime;
		public final double timeResolution;
		public final double[] slitPositions;
		public final double[] slitLengths;
		public final double[] slitWidths;

		public DetectorConfiguration(
				String cosyFile, double tiltAngleProton, double tiltAngleDeuteron, double timeResolution, double streakTime,
				double[] slitPositions, double[] slitLengths, double[] slitWidths) {
			this.cosyFile = cosyFile;
			this.tiltAngleProton = tiltAngleProton;
			this.tiltAngleDeuteron = tiltAngleDeuteron;
			this.timeResolution = timeResolution;
			this.streakTime = streakTime;
			this.slitPositions = slitPositions;
			this.slitLengths = slitLengths;
			this.slitWidths = slitWidths;
		}
	}

}
