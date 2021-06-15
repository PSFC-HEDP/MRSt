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
import util.PythonPlot;
import util.CSV;

import java.io.File;
import java.io.IOException;

/**
 * A class to handle the detector modelling
 */
public class Detector {

	private static final String SUBSTRATE_STOPPING_FILENAME = "input/stopping_power_%ss_Si.csv";
	private static final String PHOTOCATHODE_STOPPING_FILENAME = "input/stopping_power_%ss_CsI.csv";

	private static final double ESCAPE_LENGTH = 90e-4; // escape length of electrons in CsI [μm]
	private static final double ESCAPE_PROB = .70; // probability of an edge electron's escape
	private static final double WORK_FUNC = 40e-3; // work function [keV]

	private static final int INTEGRATION_RESOLUTION = 100;

	private final double substrateThickness; // [μm]
	private final double photocathodeThickness; // [μm]
	private final double bias; // [V]
	private final double meshLength; // [m]
	private final double driftLength; // [m]
	private final double dilation;
	private final double openAreaRatio;
	private final double averageGain;

	private final DiscreteFunction stoppingSi;
	private final DiscreteFunction stoppingCsI;


	/**
	 * put together the photocathode/PDDT simulacion
	 * @param ion either Particle.P or Particle.D
	 * @param substrateThickness the thickness of the silicon layer that supports
	 *                           the photocathode [μm]
	 * @param photocathodeThickness the thickness of the CsI photocathode [μm]
	 * @param bias the mean strength of the accelerating electric field [V]
	 * @param meshLength the distance across which the SE are accelerated [m]
	 * @param driftLength the length of the the PDDT [m]
	 * @param dilation the factor by which the signal is time-dilated
	 * @param openAreaRatio the fraccion of the MCT that is open
	 * @param averageGain the average number of electrons an MCT cascade produces
	 * @throws IOException if it can't find the stopping power file
	 */
	public Detector(Particle ion,
	                double substrateThickness,
	                double photocathodeThickness,
	                double bias,
	                double meshLength, double driftLength,
	                double dilation,
	                double openAreaRatio, double averageGain,
	                int integrationResolution) throws IOException {

		this.bias = bias;
		this.meshLength = meshLength;
		this.driftLength = driftLength;
		this.substrateThickness = substrateThickness;
		this.photocathodeThickness = photocathodeThickness;
		this.dilation = dilation;
		this.openAreaRatio = openAreaRatio;
		this.averageGain = averageGain;

		double[][] stoppingDataSi = CSV.read(
				new File(String.format(SUBSTRATE_STOPPING_FILENAME, ion.name)),
				',');
		stoppingSi = new DiscreteFunction(stoppingDataSi, false, true); // [keV] -> [keV/μm]
		double[][] stoppingDataCsI = CSV.read(
				new File(String.format(PHOTOCATHODE_STOPPING_FILENAME, ion.name)),
				',');
		stoppingCsI = new DiscreteFunction(stoppingDataCsI, false, true); // [keV] -> [keV/μm]
	}

	/**
	 * simulate n deuterons striking the CsI cathode and get a list of the times
	 * at which its child electrons enter the oscilliscope.
	 * @param energy the energy of the deuteron beam in MeV
	 * @return {the gain of each particle as it hits the detector,
	 *          the time at which the particle hits the detector}
	 */
	private double[][] generateDetectionEvents(int numberOfDeuterons,
	                                           double energy) throws IOException {
		double E = NumericalMethods.odeSolve( // integrate the deuteron thru the substrate
		                                      stoppingSi,
		                                      -substrateThickness,
		                                      energy*1e3,
		                                      INTEGRATION_RESOLUTION); // [keV]
		double step = photocathodeThickness/INTEGRATION_RESOLUTION; // [μm]
		double electronPerDeuteron = 0;
		for (int i = 0; i <= INTEGRATION_RESOLUTION; i ++) { // then do a fancier integral to get the number of deuterons generated in the substrate
			double x = ((double)i/INTEGRATION_RESOLUTION - 1)*photocathodeThickness; // [μm]
			double dEdx = stoppingCsI.evaluate(E);
			double dx = (i == 0 || i == INTEGRATION_RESOLUTION) ? step/2 : step;
			electronPerDeuteron += 1/WORK_FUNC*dEdx*ESCAPE_PROB*Math.exp(x/ESCAPE_LENGTH)*dx;
			E -= dEdx*step; // I'm using a first-order method here, but it should be fine because the total change in energy across the photocathode is quite small
		}

		double[][] henkeData = CSV.read(new File("input/henke_electron_data.csv"), ',');
		DiscreteFunction pdf = new DiscreteFunction(henkeData); // [eV] -> [1/eV]
		DiscreteFunction cdf = pdf.antiderivative(); // [eV] -> []
		DiscreteFunction idf = cdf.inv(); // [] -> [eV]

		final Particle e = Particle.E;
		int numberAtCathode = NumericalMethods.poisson(
				numberOfDeuterons*electronPerDeuteron);
		double[] velocityDistribution = new double[numberAtCathode]; // [J]
		for (int i = 0; i < numberAtCathode; i ++) {
			double totalElectronEnergy = idf.evaluate(Math.random()*idf.maxDatum()); // kinetick energy [eV]
			velocityDistribution[i] = Math.random()*Math.sqrt(
					2*totalElectronEnergy*(-e.charge/e.mass)); // axial velocity [m/s]
		}

		double[] timeDistribution = new double[numberAtCathode]; // [s]
		double g = -bias/meshLength*e.charge/e.mass; // [m/s^2]
		for (int i = 0; i < numberAtCathode; i ++) {
			double a = g/2;
			double b = velocityDistribution[i];
			double c = -meshLength;
			double meshTime = (-b + Math.sqrt(b*b - 4*a*c))/(2*a);
			double driftSpeed = 2*a*meshTime + b;
			timeDistribution[i] = (
					meshTime + driftLength/driftSpeed)*1e9; // [ns]
		}

		double[] gain = new double[numberAtCathode];
		for (int i = 0; i < numberAtCathode; i ++)
			gain[i] = NumericalMethods.bernoulli(openAreaRatio) ? 1 + (int)NumericalMethods.exponential(averageGain) : 0;

		return new double[][] {gain, timeDistribution};
	}

	/**
	 * @return the time it takes for an electron with no inicial energy to reach
	 * the MCT, ignoring dilation [ns]
	 */
	public double maxTime() {
		double driftSpeed = Math.sqrt(-2*bias*Particle.E.charge/Particle.E.mass); // [m/s]
		return (2*this.meshLength + this.driftLength)/driftSpeed*1e9;
	}


	public static void main(String[] args) throws IOException {
		int res = 100;
		int n = 10000;
		double[] tBounds = {-50, 5};
		Detector detector = new Detector(
				Particle.D, 100, 0.1, 1e3,
				1e-3, 1e0, 20, .60,
				1e4, 100);

		for (double energy = 11; energy < 15; energy += 1.5) {
			double[][] result = detector.generateDetectionEvents(n, energy);
			double[] gains = result[0];
			double[] times = result[1];

			double[] bins = new double[res + 1];
			double[] axis = new double[res];
			double[] dist = new double[res];
			double[] errors = new double[res];
			for (int i = 0; i < bins.length; i++)
				bins[i] = tBounds[0] + (tBounds[1] - tBounds[0])*i/res;
			for (int i = 0; i < axis.length; i++)
				axis[i] = (bins[i] + bins[i + 1])/2;
			for (int i = 0; i < times.length; i++) {
				double adjustedTime = (times[i] - detector.maxTime())/1e-3;
				int bin;
				if ((bin = NumericalMethods.bin(adjustedTime, axis)) >= 0)
					dist[bin] += gains[i];
			}
			System.out.println("Energy: " + energy);
			System.out.println("Signal amplification: " + (NumericalMethods.sum(dist)/detector.averageGain/n));
			System.out.println("Effective time resolution degradation: " + NumericalMethods.fwhm(axis, dist) + " ps");
			if (energy == 12.5)
				PythonPlot.plotLines("Energy histogram", axis, dist, errors, "Electrons");
		}
	}
}
