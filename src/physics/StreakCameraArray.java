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

import java.util.Arrays;

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
	private final double spatialResolution; // 1σ resolution (m)

	private final DiscreteFunction dxdE; // (MeV neutron -> m/MeV)
	private final DiscreteFunction dtdE; // (MeV neutron -> s/MeV)
	private final DiscreteFunction bowtieHite; // (MeV neutron -> m)

	private double[] energyBins;
	private double[] timeBins;
	private double[][] energyResponses; // (each row sums to 1)
	private double[][] timeResponses; // (each row sums to 1)

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
			 102e-6,
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
	 * @param spatialResolution the FWHM resolution of the optics/detector (m)
	 * @param optics the ion optics system (needed to set up efficient calculacion)
	 */
	public StreakCameraArray(
		  double[] slitPositions, double[] slitLengths, double[] slitWidths,
		  double sweepTime, double gain, double backgroundDensity,
		  double noiseDensity, double saturationLimitDensity,
		  double pixelArea, double spatialResolution, IonOptics optics) {
		super(backgroundDensity*Math2.gamma(4, 4, MC_RANDOM),
			  noiseDensity*Math2.gamma(4, 4, MC_RANDOM),
			  saturationLimitDensity, gain);

		if (slitPositions.length != slitWidths.length)
			throw new IllegalArgumentException("number of slits must match");
		for (double slitWidth: slitWidths)
			if (slitWidth <= 0)
				throw new IllegalArgumentException("slit widths must be positive, not " + slitWidth);

		this.slitWidths = slitWidths;

		double[] ERef = new double[FP_RESOLUTION];
		double[] xRef = new double[FP_RESOLUTION];
		double[] tRef = new double[FP_RESOLUTION];
		double[] hRef = new double[FP_RESOLUTION];
		for (int i = 0; i < FP_RESOLUTION; i ++) {
			ERef[i] = 12. + 4.*i/(FP_RESOLUTION - 1); // (MeV neutron)
			double[] rtCentral = optics.map(ERef[i]);
			xRef[i] = Math.hypot(rtCentral[0], rtCentral[2])*Math.signum(rtCentral[0]);
			tRef[i] = rtCentral[3];
			hRef[i] = 0;
			for (int k = 0; k < 1000; k ++) {
				double[] rt = optics.simulate(ERef[i], 0, false);
				if (!Double.isNaN(rt[0])) {
					if (2*Math.abs(rt[1]) > hRef[i]) hRef[i] = 2*Math.abs(rt[1]);
				}
			}
		}
		DiscreteFunction fpPosition = new DiscreteFunction(ERef, xRef, true);
		this.dxdE = fpPosition.derivative();
		DiscreteFunction fpDelay = new DiscreteFunction(ERef, tRef, true);
		this.dtdE = fpDelay.derivative();
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
		this.spatialResolution = spatialResolution/(2*Math.sqrt(2*Math.log(2))); // 1σ resolut. (m)
	}

	@Override
	public double efficiency(double energy) {
		int i = whichSlit(energy);
		if (i >= 0)
			return Math.min(1, slitWidths[i]/bowtieHite.evaluate(energy));
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
							   double[][] inSpectrum, boolean stochastic, boolean background) {
		makeSureWeHaveResponseCurves(energyBins, timeBins);

		for (double[] row : inSpectrum)
			for (double val: row)
				if (Double.isNaN(val))
					throw new IllegalArgumentException("this can't be nan.");

		double[][] outSpectrum = new double[energyBins.length-1][timeBins.length-1];
		for (int i = 0; i < energyBins.length-1; i ++) {
			double energy = (energyBins[i] + energyBins[i+1])/2;

			int slit = whichSlit(energy);
			if (slit >= 0) {
				// first construct this row according to the energy blurring
				double[] blurdSpectrum = new double[timeBins.length - 1];
				for (int k = 0; k < energyResponses[slit].length; k ++) {
					int di = k - this.energyResponses[slit].length/2;
					for (int j = 0; j < timeBins.length - 1; j ++) {
						if (i + di >= 0 && i + di < energyBins.length - 1) {
							double source = inSpectrum[i + di][j];
							double portion = energyResponses[slit][k];
							if (stochastic)
								blurdSpectrum[j] += Math2.binomial(
									  (int)source, portion, NOISE_RANDOM);
							else
								blurdSpectrum[j] += source*portion;
						}
					}
				}

				// then apply the time blurring to this row
				for (int l = 0; l < this.timeResponses[slit].length; l ++) {
					int dj = l - this.timeResponses[slit].length/2;
					for (int j = 0; j < timeBins.length - 1; j ++) {
						if (j + dj >= 0 && j + dj < timeBins.length - 1) {
							double source = blurdSpectrum[j + dj];
							double portion = this.efficiency(energy)*timeResponses[slit][l]; // accounting for quantum efficiency
							if (stochastic)
								outSpectrum[i][j] += gain*Math2.binomial(
									  (int)source, portion, NOISE_RANDOM); // as well as gain
							else
								outSpectrum[i][j] += gain*source*portion;
						}
					}
				}

				if (background) {
					for (int j = 0; j < timeBins.length - 1; j++) { // add the background
						double level = background(energy, energyBins, timeBins);
						if (stochastic) {
							double σ2 = noise(energy, energyBins, timeBins) + outSpectrum[i][j];
							outSpectrum[i][j] += Math2.normal(level, Math.sqrt(σ2), NOISE_RANDOM);
						}
						else
							outSpectrum[i][j] += level;
					}
				}
			}
		}

		return outSpectrum;
	}

	private void makeSureWeHaveResponseCurves(double[] energyBins, double[] timeBins) {
		if (this.timeResponses == null ||
			  !Arrays.equals(energyBins, this.energyBins) ||
			  !Arrays.equals(timeBins, this.timeBins)) { // (only evaluate it if it hasn't been evaluated yet)
			this.energyBins = energyBins;
			this.timeBins = timeBins;

			this.energyResponses = new double[slitWidths.length][];
			this.timeResponses = new double[slitWidths.length][];
			double energyStep = energyBins[1] - energyBins[0]; // (MeV)
			double timeStep = timeBins[1] - timeBins[0]; // (ns)

			for (int s = 0; s < this.slitWidths.length; s ++) {
				double dispersion = this.dxdE.evaluate(slitRiteBounds[s]); // (m/MeV)
				double timeSkew = this.dtdE.evaluate(slitRiteBounds[s]); // (s/MeV)

				double energyResolut = this.spatialResolution/dispersion; // (MeV)
				if (energyBins.length%2 == 1)
					energyResponses[s] = new double[energyBins.length];
				else
					energyResponses[s] = new double[energyBins.length - 1];
				for (int i = 0; i < energyResponses[s].length; i ++) {
					int di = i - energyResponses[s].length/2;
					energyResponses[s][i] = Math.exp(-Math.pow(di*energyStep/energyResolut, 2)/2.);
				}
				double energyTotal = Math2.sum(energyResponses[s]);
				for (int i = 0; i < energyResponses[s].length; i ++)
					energyResponses[s][i] /= energyTotal;

				double timeResolut = Math.hypot(
					  this.spatialResolution/streakSpeed,
					  energyResolut*timeSkew)/1e-9; // (ns)
				System.out.println("the inherent time resolution is "+this.spatialResolution/streakSpeed/1e-12+
										 " ps from the streak speed and "+energyResolut*timeSkew/1e-12+
										 " ps from the time skew");
				double[] tubeTimeResponse;
				if (timeBins.length%2 == 1)
					tubeTimeResponse = new double[timeBins.length];
				else
					tubeTimeResponse = new double[timeBins.length - 1];
				for (int j = 0; j < tubeTimeResponse.length; j ++) {
					int dj = j - tubeTimeResponse.length/2;
					tubeTimeResponse[j] = Math.exp(-Math.pow(dj*timeStep/timeResolut, 2)/2.);
				}
				double timeTotal = Math2.sum(tubeTimeResponse);
				for (int j = 0; j < tubeTimeResponse.length; j ++)
					tubeTimeResponse[j] /= timeTotal;

				double timeWidth = slitWidths[s]/streakSpeed/1e-9; // (ns)
				System.out.println("the slit degrades the resolution by "+timeWidth*1e3+"ps");
				int kernelSize = (int)Math.ceil(timeWidth/timeStep);
				if (kernelSize%2 != 1) kernelSize += 1;
				double[] slitTimeResponse = new double[kernelSize]; // bild the time response funccion kernel
				for (int j = 1; j < kernelSize - 1; j ++)
					slitTimeResponse[j] = timeStep/timeWidth;
				slitTimeResponse[0] = (1 - (kernelSize - 2)*timeStep/timeWidth)/2;
				slitTimeResponse[kernelSize - 1] = slitTimeResponse[0];
				this.timeResponses[s] = Math2.convolve(tubeTimeResponse, slitTimeResponse);
			}
		}
	}

	public double spectralRange() {
		return slitRiteBounds[slitRiteBounds.length-1] -
			  slitLeftBounds[slitLeftBounds.length-1];
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
