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
import java.util.Arrays;

import static physics.Analysis.RANDOM;

/**
 * A class to handle the detector modelling
 */
public class Detector {

	private static final String SUBSTRATE_STOPPING_FILENAME = "input/stopping_power_%ss_Si.csv";
	private static final String PHOTOCATHODE_STOPPING_FILENAME = "input/stopping_power_%ss_CsI.csv";

	private static final double ESCAPE_LENGTH = 90e-4; // escape length of electrons in CsI [μm]
	private static final double ESCAPE_PROB = .70; // probability of an edge electron's escape
	private static final double WORK_FUNC = 40e-3; // work function [keV]

	private static final double MCP_BLURRING = 0.100; // MCT response time [ns]
	private static final double CABLE_RESPONSE = 0.100; // cable response time [ns]

	private static final int EFFICIENCY_ENERGY_DEPENDENCE_RESOLUTION = 10;
	private static final int INTEGRATION_RESOLUTION = 100;
	private static final int NUM_RESPONSE_FUNCTION_TRIES = 10000;

	private final double bias; // [V]
	private final double meshLength; // [m]
	private final double driftLength; // [m]
	private final double dilation;
	private final double openAreaRatio;
	private final double averageGain;

	private final double energyFactor;

	private final DiscreteFunction electronsPerDeuteron; // [MeV] -> []

	private final DiscreteFunction henkeIDF;

	private double[] energyBins;
	private double[] timeBins;
	private double[][] responseFunction;


	/**
	 * put together the photocathode/PDDT simulacion
	 * @param ion either Particle.P or Particle.D
	 * @param substrateThickness the thickness of the silicon layer that supports
	 *                           the photocathode [μm]
	 * @param photocathodeThickness the thickness of the CsI photocathode [μm]
	 * @param photocathodeAngle the angle by which the detector plane is tilted [°]
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
	                double photocathodeThickness, double photocathodeAngle,
	                double bias,
	                double meshLength, double driftLength,
	                double dilation,
	                double openAreaRatio, double averageGain,
	                int integrationResolution) throws IOException {

		double A = ion.mass/Particle.N.mass;
		this.energyFactor = 4*A/Math.pow(A + 1, 2);

		this.bias = bias;
		this.meshLength = meshLength;
		this.driftLength = driftLength;
		this.dilation = dilation;
		this.openAreaRatio = openAreaRatio;
		this.averageGain = averageGain;

		double[][] stoppingDataSi = CSV.read(
				new File(String.format(SUBSTRATE_STOPPING_FILENAME, ion.name)),
				',');
		DiscreteFunction stoppingSi = new DiscreteFunction(stoppingDataSi, false, true); // [keV] -> [keV/μm]
		double[][] stoppingDataCsI = CSV.read(
				new File(String.format(PHOTOCATHODE_STOPPING_FILENAME, ion.name)),
				',');
		DiscreteFunction stoppingCsI = new DiscreteFunction(stoppingDataCsI, false, true); // [keV] -> [keV/μm]

		double[] deuteronEnergies = new double[EFFICIENCY_ENERGY_DEPENDENCE_RESOLUTION+1];
		double[] electronsPerDeuteron = new double[EFFICIENCY_ENERGY_DEPENDENCE_RESOLUTION+1];
		for (int i = 0; i < deuteronEnergies.length; i ++) {
			deuteronEnergies[i] = 9. + (15. - 9.)*i/deuteronEnergies.length;
			double E = NumericalMethods.odeSolve( // integrate the deuteron thru the substrate
					stoppingSi,
					- substrateThickness,
					deuteronEnergies[i]*1e3,
					INTEGRATION_RESOLUTION); // [keV]
			double step = photocathodeThickness/INTEGRATION_RESOLUTION; // [μm]
			for (int j = 0; j <= INTEGRATION_RESOLUTION; j ++) { // then do a fancier integral to get the number of deuterons generated in the substrate
				double x = ((double) j/INTEGRATION_RESOLUTION - 1)*photocathodeThickness; // [μm]
				double dEdx = stoppingCsI.evaluate(E)/Math.cos(Math.toRadians(photocathodeAngle));
				double dx = (j == 0 || j == INTEGRATION_RESOLUTION) ? step/2 : step;
				electronsPerDeuteron[i] += 1/WORK_FUNC*dEdx*ESCAPE_PROB*Math.exp(x/ESCAPE_LENGTH)*dx;
				E -= dEdx*step; // I'm using a first-order method here, but it should be fine because the total change in energy across the photocathode is quite small
			}
		}
		this.electronsPerDeuteron = new DiscreteFunction(deuteronEnergies,
		                                                 electronsPerDeuteron);

		double[][] henkeData = CSV.read(new File("input/henke_electron_data.csv"), ',');
		DiscreteFunction henkePDF = new DiscreteFunction(henkeData); // [eV] -> [1/eV]
		DiscreteFunction henkeCDF = henkePDF.antiderivative(); // [eV] -> []
		this.henkeIDF = henkeCDF.inv(); // [] -> [eV]
	}


	/**
	 * compute the detected spectrum given a deuteron spectrum at the photocathode
	 * @param stochastic whether to apply noise to the result
	 */
	public double[][] response(double[] energyBins, double[] timeBins,
	                           double[][] inSpectrum, boolean stochastic) {
		if (this.responseFunction == null ||
				!Arrays.equals(energyBins, this.energyBins) ||
				!Arrays.equals(timeBins, this.timeBins)) // the full nmxnm believed transfer matrix
			evaluateResponseFunction(energyBins, timeBins);

		double[][] outSpectrum = new double[energyBins.length-1][timeBins.length-1];
		for (int i = 0; i < energyBins.length-1; i ++)
			for (int j = 0; j < timeBins.length-1; j ++)
				for (int k = 0; k < timeBins.length-1; k ++)
					outSpectrum[i][j] += inSpectrum[i][k]*responseFunction[i][j - k + responseFunction[0].length/2];

		if (stochastic) {
			for (int i = 0; i < energyBins.length-1; i ++)
				for (int j = 0; j < timeBins.length-1; j ++)
					outSpectrum[i][j] = NumericalMethods.erlang(
							NumericalMethods.poisson(outSpectrum[i][j], RANDOM),
							averageGain, RANDOM);
		}

		return outSpectrum;
	}


	/**
	 * the number of signal electrons created for every incident deuteron
	 * @param energy energy of the deuteron [MeV]
	 */
	public double efficiency(double energy) {
		return this.electronsPerDeuteron.evaluate(energy)*this.openAreaRatio;
	}


	/**
	 * the average number of cable electrons from each signal electron
	 */
	public double gain() {
		return this.averageGain;
	}


	/**
	 * evaluate the response funccion of the detector and save it.
	 * @param energyBins the binning by which to encode the energy dependence [MeV]
	 * @param timeBins the binning by which to encode the response funccion [ns]
	 */
	private void evaluateResponseFunction(double[] energyBins, double[] timeBins) {
		int n = energyBins.length - 1;
		int m = timeBins.length - 1;

		double[] responseTimeBins = new double[2*m]; // create a new time axis that is twice as big
		double[] pddtResponse = new double[2*m - 1]; // and compute the temporal response, ignoring energy

		final Particle e = Particle.E;
		double expectedTime = (2*meshLength + driftLength)/Math.sqrt(2*bias*(-e.charge/e.mass))/1e-9;
		for (int j = 0; j < 2*m; j ++)
			responseTimeBins[j] = expectedTime + (j - m + 0.5)*(timeBins[1] - timeBins[0]);

		for (int k = 0; k < NUM_RESPONSE_FUNCTION_TRIES; k ++) {
			double initialTime = (RANDOM.nextDouble() - 0.5)*(timeBins[1] - timeBins[0]);
			double electronEnergy = henkeIDF.evaluate(Math.random()*henkeIDF.maxDatum());
			double a = bias/meshLength*(-e.charge/e.mass)/2; // [m/s^2]
			double b = RANDOM.nextDouble()*Math.sqrt(
					2*electronEnergy*(-e.charge/e.mass)); // [m/s]
			double c = -meshLength; // [m]
			double meshTime = (-b + Math.sqrt(b*b - 4*a*c))/(2*a);
			double driftSpeed = 2*a*meshTime + b;
			double totalTime = initialTime + (meshTime + driftLength/driftSpeed)*1e9; // [ns]
			int j = NumericalMethods.bin(totalTime, responseTimeBins);
			if (j >= 0)
				pddtResponse[j] += openAreaRatio/NUM_RESPONSE_FUNCTION_TRIES;
		}

		double[] cableResponse = new double[2*m - 1]; // TODO use actual data
		double sum = 0;
		for (int j = 0; j < cableResponse.length; j ++) {
			double t = (j + 0.5)*(timeBins[1] - timeBins[0]);
			cableResponse[j] = t*Math.exp(-t/5e-3); // use a gamma distribution
			sum += cableResponse[j];
		}
		for (int j = 0; j < cableResponse.length; j ++)
			cableResponse[j] /= sum;

		double[] monoenergeticResponse = new double[2*m - 1]; // convolve the two components
		for (int j = 0; j < 2*m - 1; j ++)
			for (int k = j; k < 2*m - 1; k ++)
				monoenergeticResponse[k] += cableResponse[k - j] * pddtResponse[j];

		this.responseFunction = new double[n][2*m - 1]; // combine with the energy dependence of the efficiency
		for (int i = 0; i < n; i ++) {
			double energy = (energyBins[i] + energyBins[i + 1])/2.;
			double photocathodeEfficiency = electronsPerDeuteron.evaluate(energy*energyFactor);
			for (int j = 0; j < 2*m - 1; j ++)
				responseFunction[i][j] = photocathodeEfficiency*monoenergeticResponse[j];
		}

		this.energyBins = energyBins; // save these energy and time bins
		this.timeBins = timeBins;
	}


	/**
	 * compute and return the response function of the MCP and faraday cups
	 * @param times the time axis for the funccion, starting at 0 [ns]
	 * @return the voltage at each time after a delta funccion of signal electrons
	 */
	private double[] MCPResponseFunction(double[] times) {
		double[] intermediate = new double[times.length];
		for (int i = 0; i < times.length; i ++)
			intermediate[i] = Math.exp(-times[i]/(CABLE_RESPONSE/dilation));//1/(CABLE_RESPONSE/dilation + times[i]);
		double[] convolved = new double[times.length];
		for (int i = 0; i < times.length; i ++)
			for (int j = 0; j < i; j ++)
				convolved[i] += Math.exp(-(times[i] - times[j])/(MCP_BLURRING/dilation))*intermediate[j]/1e8;
//				convolved[i] += Math.exp(-Math.pow(times[i] - times[j] - 3*MCP_BLURRING/dilation, 2)/
//						                         (2*Math.pow(MCP_BLURRING/dilation, 2)))*intermediate[j]/1e11;
		return convolved;
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
		double[] tBounds = {-60, 60};
		Detector detector = new Detector(
				Particle.D, 100, 0.1, 68, 1e3,
				1e-3, 1e0, 20, .60,
				1e4, 100);

		for (double energy = 11; energy < 15; energy += 1.5) {
			double[] bins = new double[res + 1];
			double[] axis = new double[res];
			double[] errors = new double[res];
			for (int i = 0; i < bins.length; i ++)
				bins[i] = tBounds[0] + (tBounds[1] - tBounds[0])*i/res;
			for (int i = 0; i < axis.length; i ++)
				axis[i] = (bins[i] + bins[i + 1])/2;

			detector.evaluateResponseFunction(new double[] {10, 14}, bins);
			double[] dist = detector.responseFunction[0];

			double[] kernelAxis = new double[res];
			for (int i = 0; i < kernelAxis.length; i ++)
				kernelAxis[i] = (axis[i] - axis[0])/1e3;
			double[] kernel = detector.MCPResponseFunction(kernelAxis);
			double[] blurredDist = new double[res];
			for (int i = 0; i < res; i ++)
				for (int j = 0; j <= i; j ++)
					blurredDist[i] += kernel[j]*dist[i - j];

			System.out.println("Energy: " + energy);
			System.out.println("Signal amplification: " + (NumericalMethods.sum(dist)/detector.averageGain/n));
			System.out.println("Effective time resolution degradation: " + NumericalMethods.fwhm(axis, dist) + " ps");
			if (energy == 12.5) {
//				PythonPlot.plotLines("Energy histogram", axis, "Time [ps]", dist, errors, "Electrons");
//				PythonPlot.plotLines("MCP response function", kernelAxis, "Time [ps]", kernel, errors, " ");
				PythonPlot.plotLines("Measured signal", axis, "Time [ps]", blurredDist, errors, "Signal level");
			}
		}
	}
}
