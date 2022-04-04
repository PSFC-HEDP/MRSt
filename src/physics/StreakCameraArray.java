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
import util.Math2.DiscreteFunction;

import static physics.Analysis.MC_RANDOM;
import static physics.Analysis.NOISE_RANDOM;


public class StreakCameraArray extends Detector {

	private static final int FP_RESOLUTION = 50;

	private final double[] slitWidths; // (m)
	private final double streakSpeed; // (m/s in photocathode dimensions)
	private final double[] slitLeftBounds; // (MeV neutron)
	private final double[] slitRiteBounds; // (MeV neutron)
	private final double gain; // (counts/deuteron)
	private final double pixelArea; // (m^2)

	private final DiscreteFunction dxdE; // (MeV neutron -> m/MeV)
	private final DiscreteFunction bowtieHite; // (MeV neutron -> m)

	public StreakCameraArray(
		  DetectorConfiguration config, IonOptics ionOptis) {
		this(config.slitPositions,
			 config.slitLengths,
			 config.slitWidths,
			 config.streakTime,
			 2.4/Math.cos(Math.toRadians(config.tiltAngle)) * 51,
			 81,
			 81,
			 40_000,
			 25e-6*25e-6,
			 ionOptis);
	}

	/**
	 * create a new streak camera array to be used with the ion optics
	 * @param slitLengths the length of each slit in the dispersion direccion (m)
	 * @param slitWidths the height of each slit in the nondispersion direccion (m)
	 * @param sweepTime the amount of time it takes to sweep (s)
	 * @param backgroundDensity the inherent background level in the detector (counts/pixel)
	 * @param noiseDensity the inherent noise level in the detector (counts^2/pixel)
	 * @param slitPositions the centers of the slits on the focal plane
	 * @param optics the ion optics system (needed to set up efficient calculacion)
	 */
	public StreakCameraArray(
		  double[] slitPositions, double[] slitLengths, double[] slitWidths,
		  double sweepTime, double gain, double backgroundDensity,
		  double noiseDensity, double saturationLimitDensity,
		  double pixelArea, IonOptics optics) {
		super(backgroundDensity*Math2.gamma(4, 4, MC_RANDOM),
			  noiseDensity*Math2.gamma(4, 4, MC_RANDOM),
			  saturationLimitDensity);

		if (slitPositions.length != slitWidths.length)
			throw new IllegalArgumentException("number of slits must match");
		for (double slitWidth: slitWidths)
			if (slitWidth <= 0)
				throw new IllegalArgumentException("slit widths must be positive, not " + slitWidth);

		this.slitWidths = slitWidths;

		double[] ERef = new double[FP_RESOLUTION];
		double[] xRef = new double[FP_RESOLUTION];
		double[] hRef = new double[FP_RESOLUTION];
		for (int i = 0; i < FP_RESOLUTION; i ++) {
			ERef[i] = 12. + 4.*i/(FP_RESOLUTION - 1); // (MeV neutron)
			double num = 0, xSum = 0, yMax = 0;
			for (int k = 0; k < 1000; k ++) {
				double[] rt = optics.map(ERef[i]);
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
		this.dxdE = fpPosition.derivative();
		DiscreteFunction fpEnergy = fpPosition.inv().indexed(FP_RESOLUTION);
		this.bowtieHite = new DiscreteFunction(ERef, hRef, true);

		this.slitLeftBounds = new double[slitPositions.length];
		this.slitRiteBounds = new double[slitPositions.length];
		for (int i = 0; i < this.slitLeftBounds.length; i ++) {
			this.slitLeftBounds[i] = fpEnergy.evaluate(
				  slitPositions[i] - slitLengths[i]/2);
			this.slitRiteBounds[i] = fpEnergy.evaluate(
				  slitPositions[i] + slitLengths[i]/2);
		}

		this.streakSpeed = Math2.min(slitLengths)/sweepTime; // (m/s in photocathode scale)
		this.gain = gain;
		this.pixelArea = pixelArea;
	}

	@Override
	public double gain(double energy) {
		int i = whichSlit(energy);
		if (i >= 0)
			return this.gain*Math.min(1, slitWidths[i]/bowtieHite.evaluate(energy));
		return 0;
	}

	@Override
	public double pixelsPerBin(double energy, double[] energyBins, double[] timeBins) {
		double binSize = (energyBins[1] - energyBins[0]) * (timeBins[1] - timeBins[0]); // (MeV*ns/bin)
		double dispersion = dxdE.evaluate(energy)*streakSpeed*1e-9; // (m^2/(MeV*ns))
		return binSize*dispersion/pixelArea;
	}

	@Override
	public double[][] response(double[] energyBins, double[] timeBins,
							   double[][] inSpectrum, boolean stochastic) {
		for (double[] row : inSpectrum)
			for (double val: row)
				if (Double.isNaN(val))
					throw new IllegalArgumentException("this can't be nan.");

		double[][] timeResponses = new double[slitWidths.length][];
		for (int i = 0; i < slitWidths.length; i ++) {
			double timeWidth = slitWidths[i]/streakSpeed/1e-9; // (ns)
			double timeStep = timeBins[1] - timeBins[0]; // (ns)
			int kernelSize = (int) Math.ceil(timeWidth/timeStep);
			if (kernelSize%2 != 1)
				kernelSize += 1;
			timeResponses[i] = new double[kernelSize]; // bild the time response funccion kernel
			for (int j = 1; j < kernelSize - 1; j ++)
				timeResponses[i][j] = timeStep/timeWidth;
			timeResponses[i][0] = (1 - (kernelSize - 2)*timeStep/timeWidth)/2;
			timeResponses[i][kernelSize - 1] = timeResponses[i][0];
		}

		double[][] outSpectrum = new double[energyBins.length-1][timeBins.length-1];
		for (int i = 0; i < energyBins.length-1; i ++) {
			double energy = (energyBins[i] + energyBins[i+1])/2;
			int slit = whichSlit(energy);
			if (slit >= 0) {
				for (int j = 0; j < timeBins.length - 1; j ++) // add the background
					outSpectrum[i][j] = background(energy, energyBins, timeBins);

				double gain = this.gain(energy); // then convolve in the signal
				for (int j = 0; j < timeBins.length - 1; j ++) {
					for (int l = 0; l < timeResponses[slit].length; l ++) {
						int dj = l - timeResponses[slit].length/2;
						if (j + dj >= 0 && j + dj < timeBins.length - 1)
							outSpectrum[i][j] += gain*timeResponses[slit][l]*inSpectrum[i][j + dj];
					}
				}
				if (stochastic) { // add noise
					for (int j = 0; j < timeBins.length - 1; j ++) {
						double σ = Math.sqrt(noise(energy, energyBins, timeBins));
						if (outSpectrum[i][j] > 0)
							outSpectrum[i][j] = Math2.normal(outSpectrum[i][j], σ, NOISE_RANDOM);
					}
				}
			}
		}

		return outSpectrum;
	}

	/**
	 * find the slit that corresponds to this energy
	 * @param energy the neutron energy in MeV
	 * @return the index of the slit that gets it, or -1 if it is in no slit
	 */
	private int whichSlit(double energy) {
		for (int i = 0; i < slitLeftBounds.length; i ++)
			if (energy >= slitLeftBounds[i] && energy <= slitRiteBounds[i])
				return i;
		return -1;
	}
}
