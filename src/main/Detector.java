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
package main;

import main.NumericalMethods.DiscreteFunction;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;

/**
 * A class to handle the detector modelling
 */
public class Detector {

	public static final double BIAS = 1e3; // [V]
	public static final double MESH_LENGTH = .001; // [m]
	public static final double DRIFT_LENGTH = 1; // [m]
	public static final double SUBSTRATE_THICKNESS = 100; // [μm]
	public static final double PHOTOCATHODE_THICKNESS = 0.1; // [μm]
	public static final double DILATION = 20;
	public static final double OPEN_AREA_RATIO = .60;
	public static final double AVERAGE_GAIN = 1e4;
	private static final int RESOLUTION = 36;
	private static final int INTEGRATION_RESOLUTION = 100;
	private static final int DISTRIBUTION_RESOLUTION = 100;

	/**
	 * simulate n deuterons striking the CsI cathode and get a list of the times
	 * at which its child electrons enter the oscilliscope.
	 * @param energy the energy of the deuteron beam in MeV
	 * @return {the gain of each particle as it hits the detector,
	 *          the time at which the particle hits the detector}
	 */
	private static double[][] generateDetectionEvents(int numberOfDeuterons,
	                                                  double energy) throws IOException {
		double[][] stoppingSi;
		double[][] stoppingCsI;
		try {
			stoppingSi = CSV.read(new File("input/stopping_power_deuterons_Si.csv"), ',');
			stoppingCsI = CSV.read(new File("input/stopping_power_deuterons_CsI.csv"), ',');
		} catch (IOException e) {
			stoppingSi = new double[8][2];
			stoppingCsI = new double[8][2];
			e.printStackTrace();
		}
		DiscreteFunction dEdxSi = new DiscreteFunction(stoppingSi, false, true); // [keV] -> [keV/μm]
		DiscreteFunction dEdxCsI = new DiscreteFunction(stoppingCsI, false, true); // [keV] -> [keV/μm]
		double E = NumericalMethods.odeSolve( // integrate the deuteron thru the substrate
				dEdxSi,
				-SUBSTRATE_THICKNESS,
				energy*1e3,
				INTEGRATION_RESOLUTION); // [keV]
		double step = PHOTOCATHODE_THICKNESS/INTEGRATION_RESOLUTION; // [μm]
		double P0 = .70;
		double ɛ = 40e-3; // work function [keV]
		double L = 90e-4; // [μm]
		double electronPerDeuteron = 0;
		for (int i = 0; i <= INTEGRATION_RESOLUTION; i ++) { // then do a fancier integral to get the number of deuterons generated in the substrate
			double x = ((double)i/INTEGRATION_RESOLUTION - 1)*PHOTOCATHODE_THICKNESS; // [μm]
			double dEdx = dEdxCsI.evaluate(E);
			double dx = (i == 0 || i == INTEGRATION_RESOLUTION) ? step/2 : step;
			electronPerDeuteron += 1/ɛ*dEdx*P0*Math.exp(x/L)*dx;
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
			double totalElectronEnergy = idf.evaluate(Math.random()*idf.maxDatum()); // [eV]
			double axialElectronEnergy = Math.random()*totalElectronEnergy; // [eV]
			velocityDistribution[i] = Math.sqrt(
					2*axialElectronEnergy*(-e.charge/e.mass)); // [m/s]
		}

		double[] timeDistribution = new double[numberAtCathode]; // [s]
		double g = -BIAS/MESH_LENGTH*e.charge/e.mass; // [m/s^2]
		for (int i = 0; i < numberAtCathode; i ++) {
			double a = g/2;
			double b = velocityDistribution[i];
			double c = -MESH_LENGTH;
			timeDistribution[i] = (-b + Math.sqrt(b*b - 4*a*c))/(2*a);
			double vDrift = 2*a*timeDistribution[i] + b;
			timeDistribution[i] += DRIFT_LENGTH/vDrift;
			timeDistribution[i] *= DILATION*1e9; // dilate and convert to ns
		}

		double[] gain = new double[numberAtCathode];
		for (int i = 0; i < numberAtCathode; i ++)
			gain[i] = NumericalMethods.bernoulli(OPEN_AREA_RATIO) ? 1 + (int)NumericalMethods.exponential(AVERAGE_GAIN) : 0;

		return new double[][] {gain, timeDistribution};
	}

	/**
	 * send 1D data to a Python script for plotting in MatPlotLib
	 * @param x the data for the x axis
	 * @param yDatums {data for the y axis, error bar width for the y axis, label for that data,
	 *   ...}
	 * @throws IOException if there's an issue talking to disk
	 */
	private static void plotLines(String name, double[] x, Object... yDatums) throws IOException, InterruptedException {
		double[][] ys = new double[yDatums.length/3][];
		double[][] Δs = new double[yDatums.length/3][];
		String[] yLabels = new String[yDatums.length/3];
		for (int i = 0; i < yDatums.length/3; i ++) {
			ys[i] = (double[]) yDatums[3*i];
			Δs[i] = (double[]) yDatums[3*i + 1];
			yLabels[i] = (String) yDatums[3*i + 2];
		}

		new File("output/").mkdir();
		CSV.writeColumn(x, new File(String.format("output/%s_x.csv", "data")));
		for (int i = 0; i < ys.length; i ++) {
			CSV.writeColumn(ys[i], new File(String.format("output/%s_y_%d.csv", "Data", i)));
			CSV.writeColumn(Δs[i], new File(String.format("output/%s_err_%d.csv", "Data", i)));
		}
		ProcessBuilder plotPB = new ProcessBuilder("python", "src/python/plot1.py",
		                                           "Time (ns)", String.join("\n", yLabels),
		                                           "data", name, Integer.toString(ys.length));
		Process p = plotPB.start();
		System.out.println("finishd plotting with exit code "+p.waitFor());
	}


	public static void main(String[] args) throws IOException, InterruptedException {
		int n = 10000;
		double[][] result = generateDetectionEvents(n, 12.5);
		double[] gains = result[0];
		double[] times = result[1];

		double[] bins = new double[RESOLUTION + 1];
		double[] axis = new double[RESOLUTION];
		double[] dist = new double[RESOLUTION];
		double[] errors = new double[RESOLUTION];
		for (int i = 0; i < bins.length; i ++)
			bins[i] = 1065 + 4.*i/RESOLUTION;
		for (int i = 0; i < axis.length; i ++)
			axis[i] = (bins[i] + bins[i+1])/2;
		for (int i = 0; i < times.length; i ++) {
			int bin;
			if ((bin = NumericalMethods.bin(times[i], axis)) >= 0)
				dist[bin] += gains[i];
		}
		System.out.println("Signal amplification: "+(NumericalMethods.sum(dist)/AVERAGE_GAIN/n));
		System.out.println("Effective time resolution degradation: "+(NumericalMethods.fwhm(axis, dist)/DILATION*1e3)+" ps");
		plotLines("Energy histogram", axis, dist, errors, "Electrons");
	}
}
