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

import util.NumericalMethods;
import util.NumericalMethods.DiscreteFunction;

import static physics.Analysis.NOISE_RANDOM;


public class StreakCameraArray implements Detector {

	private static final int FP_RESOLUTION = 30;

	private final double slitWidth; // (m)
	private final double streakSpeed; // (m/s in photocathode dimensions)
	private final double[] slitLeftBounds; // (MeV neutron)
	private final double[] slitRiteBounds; // (MeV neutron)
	private final double gain; // (counts/deuteron)
	private final double noiseDensity; // (counts^2/m^2)
	private final double backgroundDensity; // (counts/m^2)

	private final DiscreteFunction dxdE; // (MeV neutron -> m/MeV)
	private final DiscreteFunction bowtieHite; // (MeV neutron -> m)

	/**
	 * create a new streak camera array to be used with the ion optics
	 * @param slitLength the length of each slit in the dispersion direccion (m)
	 * @param slitWidth the height of each slit in the nondispersion direccion (m)
	 * @param sweepTime the amount of time it takes to sweep (s)
	 * @param noiseDensity the inherent noise level in the detector (counts^2/m^2)
	 * @param backgroundDensity the inherent background level in the detector (counts/m^2)
	 * @param slitPositions the centers of the slits on the focal plane
	 * @param optics the ion optics system (needed to set up efficient calculacion)
	 */
	public StreakCameraArray(
		  double slitLength, double slitWidth,
		  double sweepTime, double gain, double noiseDensity, double backgroundDensity,
		  double[] slitPositions,
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
		this.dxdE = fpPosition.derivative();
		DiscreteFunction fpEnergy = fpPosition.inv().indexed(FP_RESOLUTION);
		this.bowtieHite = new DiscreteFunction(ERef, hRef, true);

		double fpCenter = fpPosition.evaluate(14);
		for (int i = 0; i < slitPositions.length; i ++)
			slitPositions[i] -= fpCenter;
		this.slitLeftBounds = new double[slitPositions.length];
		this.slitRiteBounds = new double[slitPositions.length];
		for (int i = 0; i < this.slitLeftBounds.length; i ++) {
			this.slitLeftBounds[i] = fpEnergy.evaluate(
				  slitPositions[i] - slitLength/2);
			this.slitRiteBounds[i] = fpEnergy.evaluate(
				  slitPositions[i] + slitLength/2);
		}

		this.streakSpeed = slitLength/sweepTime; // (m/s in photocathode scale)
		this.gain = gain;
		this.backgroundDensity = backgroundDensity;//*NumericalMethods.gamma(4, 4, MC_RANDOM); // (counts/ns*MeV)
		this.noiseDensity = noiseDensity;
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
		return this.gain;
	}

	@Override
	public double noise(double energy, double[] energyBins, double[] timeBins) {
		double binScale = dxdE.evaluate(energy)*streakSpeed*1e-9; // (m^2/(MeV*ns))
		return noiseDensity*binScale
			  *(energyBins[1] - energyBins[0])*(timeBins[1] - timeBins[0]);
	}

	@Override
	public double background(double energy, double[] energyBins, double[] timeBins) {
		double binScale = dxdE.evaluate(energy)*streakSpeed*1e-9; // (m^2/(MeV*ns))
		return backgroundDensity*binScale
			  *(energyBins[1] - energyBins[0])*(timeBins[1] - timeBins[0]);
	}

	@Override
	public double[][] response(double[] energyBins, double[] timeBins,
							   double[][] inSpectrum, boolean stochastic) {
		double timeWidth = slitWidth/streakSpeed/1e-9; // (ns)
		double timeStep = timeBins[1] - timeBins[0]; // (ns)
		int kernelSize = 2 + (int)Math.ceil(timeWidth/timeStep);
		if (kernelSize%2 != 1)
			kernelSize += 1;
		double[] timeResponse = new double[kernelSize]; // bild the time response funccion kernel
		for (int i = 1; i < kernelSize - 1; i ++)
			timeResponse[i] = timeStep/timeWidth;
		timeResponse[0] = timeResponse[kernelSize-1] = (timeWidth - (kernelSize - 2)*timeStep)/timeWidth/2;

		double[][] outSpectrum = new double[energyBins.length-1][timeBins.length-1];
		for (int i = 0; i < energyBins.length-1; i ++) {
			double energy = (energyBins[i] + energyBins[i+1])/2;
			double efficiency = this.efficiency(energy);
			if (efficiency != 0) {
				for (int j = 0; j < timeBins.length - 1; j ++)
					outSpectrum[i][j] = background(energy, energyBins, timeBins);
				for (int j = 0; j < timeBins.length - 1; j++) {
					for (int l = 0; l < kernelSize; l++) {
						int dj = l - kernelSize/2;
						if (j + dj >= 0 && j + dj < timeBins.length - 1)
							outSpectrum[i][j] += efficiency*gain*timeResponse[l]*inSpectrum[i][j + dj];
					}
				}
				if (stochastic) {
					for (int j = 0; j < timeBins.length - 1; j ++) {
						double σ = Math.sqrt(
							  noise(energy, energyBins, timeBins) +
							  (outSpectrum[i][j] - background(energy, energyBins, timeBins))/gain);
						if (outSpectrum[i][j] > 0)
							outSpectrum[i][j] = NumericalMethods.normal(outSpectrum[i][j], σ, NOISE_RANDOM);
					}
				}
			}
		}

		return outSpectrum;
	}
}
