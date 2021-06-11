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

import java.io.File;
import java.io.IOException;
import java.util.Arrays;

/**
 * A class to handle the detector modelling
 */
public class Detector {

	public static double BIAS = 2e3; // [V]
	public static double LENGTH = 1; // [m]
	public static double THICKNESS = 1e3; // [m]
	public static double DILATION = 20;
	public static double OPEN_AREA_RATIO = 60e-2;
	public static double AVERAGE_GAIN = 10;
	private static final int RESOLUTION = 36;
	private static final int INTEGRATION_RESOLUTION = 100;

	/**
	 * simulate n deuterons striking the CsI cathode and get a list of the times
	 * at which its child electrons enter the oscilliscope.
	 * @param energy the energy of the deuteron beam in eV
	 * @return {the gain of each particle as it hits the detector,
	 *          the time at which the particle hits the detector}
	 */
	private static double[][] generateDetectionEvents(int n, double energy) {
		double[] intEnergy = new double[INTEGRATION_RESOLUTION]; // TODO: get stopping power data for CsI

		double electronPerDeuteron = 1.8;
		final Particle e = Particle.E;
		int numberAtCathode = (int)Math.round(NumericalMethods.poisson(
				n*electronPerDeuteron)); // TODO: calculate energy-dependent mean
		double[] energyDistribution = new double[numberAtCathode]; // [J]
		for (int i = 0; i < numberAtCathode; i ++)
			energyDistribution[i] = NumericalMethods.exponential(-0.5*e.charge) + NumericalMethods.exponential(-0.5*e.charge); // TODO: get the actual distribution from Henke (1981)

		double[] timeDistribution = new double[numberAtCathode]; // [s]
		double g = -BIAS/LENGTH*e.charge/e.mass; // [m/s^2]
		for (int i = 0; i < numberAtCathode; i ++) {
			double a = g/2;
			double b = Math.sqrt(2*energyDistribution[i]/e.mass); // TODO: is this energy distribucion monodireccional?
			double c = -LENGTH;
			timeDistribution[i] = DILATION*1e9*(-b + Math.sqrt(b*b - 4*a*c))/(2*a);
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
		System.out.println(p.waitFor());
	}


	public static void main(String[] args) throws IOException, InterruptedException {
		double[][] result = generateDetectionEvents(10000, 12.5);
		double[] gains = result[0];
		double[] times = result[1];

		double[] bins = new double[RESOLUTION + 1];
		double[] axis = new double[RESOLUTION];
		double[] dist = new double[RESOLUTION];
		double[] errors = new double[RESOLUTION];
		for (int i = 0; i < bins.length; i ++)
			bins[i] = 1420 + 100.*i/RESOLUTION;
		for (int i = 0; i < axis.length; i ++)
			axis[i] = (bins[i] + bins[i+1])/2;
		for (int i = 0; i < times.length; i ++) {
			int bin;
			if ((bin = NumericalMethods.bin(times[i], axis)) >= 0)
				dist[bin] += gains[i];
		}
		plotLines("Energy histogram", axis, dist, errors, "Electrons");
	}
}
