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
package util;

import java.io.File;
import java.io.IOException;

/**
 * a utility class for using my python plotting scripts
 */
public class PythonPlot {
	/**
	 * send 2D data to a Python script for plotting in MatPlotLib
	 * @throws IOException if there's an issue talking to disc
	 */
	public static void plotHeatmap(
			double[] x, double[] y, double[][] z,
			String xLabel, String yLabel, String title) throws IOException {
		new File("output/").mkdir();
		CSV.writeColumn(x, new File(String.format("output/%s_x.csv", title)));
		CSV.writeColumn(y, new File(String.format("output/%s_y.csv", title)));
		CSV.write(z, new File(String.format("output/%s_z.csv", title)), ',');
		ProcessBuilder plotPB = new ProcessBuilder("python", "src/python/plot2.py",
		                                           xLabel, yLabel, title);
		plotPB.start();
	}


	/**
	 * send 1D data to a Python script for plotting in MatPlotLib
	 * @param x the data for the x axis
	 * @param yDatums {data for the y axis, error bar width for the y axis, label
	 *                for that data, ...}
	 * @throws IOException if there's an issue talking to disk
	 */
	public static void plotLines(String title, String name, double[] x, String xLabel,
								 Object... yDatums) throws IOException {
		double[][] ys = new double[yDatums.length/3][];
		double[][] Δs = new double[yDatums.length/3][];
		String[] yLabels = new String[yDatums.length/3];
		for (int i = 0; i < yDatums.length/3; i ++) {
			ys[i] = (double[]) yDatums[3*i];
			Δs[i] = (double[]) yDatums[3*i+1];
			yLabels[i] = (String) yDatums[3*i+2];
			if (Δs[i] == null)
				Δs[i] = new double[ys[i].length];
		}

		new File("output/").mkdir();
		CSV.writeColumn(x, new File(String.format("output/%s_x.csv", title)));
		for (int i = 0; i < ys.length; i ++) {
			CSV.writeColumn(ys[i], new File(String.format("output/%s_y_%d.csv", title, i)));
			CSV.writeColumn(Δs[i], new File(String.format("output/%s_err_%d.csv", title, i)));
		}
		ProcessBuilder plotPB = new ProcessBuilder("python", "src/python/plot1.py",
		                                           xLabel, String.join("\n", yLabels),
		                                           title, name, Integer.toString(ys.length));
		plotPB.start();
	}


	/**
	 * send 1D data to a Python script for plotting in MatPlotLib
	 * @throws IOException if there's an issue talking to disk
	 */
	public static void compareHeatmap(double[] x, double[] y, double[][] z0, double[][] z1,
	                                  String xLabel, String yLabel, String title0, String title1) throws IOException {
		new File("output/").mkdir();
		CSV.writeColumn(x, new File(String.format("output/%s_x.csv", title0)));
		CSV.writeColumn(y, new File(String.format("output/%s_y.csv", title0)));
		CSV.write(z0, new File(String.format("output/%s_z.csv", title0)), ',');
		CSV.write(z1, new File(String.format("output/%s_z.csv", title1)), ',');
		ProcessBuilder plotPB = new ProcessBuilder("python", "src/python/compare2.py",
		                                           xLabel, yLabel, title0, title1);
		plotPB.start();
	}


	/**
	 * send 1D data to a Python script for plotting in MatPlotLib
	 * @throws IOException if there's an issue talking to disk
	 */
	public static void compareHeatmap(double[] x, double[] y,
	                                  double[][] z0, double[][] z1, double[][] z2,
	                                  String xLabel, String yLabel,
	                                  String title0, String title1, String title2) throws IOException {
		new File("output/").mkdir();
		CSV.writeColumn(x, new File(String.format("output/%s_x.csv", title0)));
		CSV.writeColumn(y, new File(String.format("output/%s_y.csv", title0)));
		CSV.write(z0, new File(String.format("output/%s_z.csv", title0)), ',');
		CSV.write(z1, new File(String.format("output/%s_z.csv", title1)), ',');
		CSV.write(z2, new File(String.format("output/%s_z.csv", title2)), ',');
		ProcessBuilder plotPB = new ProcessBuilder("python", "src/python/compare2.py",
		                                           xLabel, yLabel, title0, title1, title2);
		plotPB.start();
	}


	/**
	 * show a graphic of the focal plane given these measurements
	 * @param energies the energies at which particles were generated
	 * @param particles for each energy, a flattend array of random {x, y, t}
	 * @param spectrum the spectral density in each energy bin
	 * @param slitPositions the position of the center of each slit (m)
	 * @param slitLengths the length of each slit (m)
	 * @param slitWidths the width of each slit (m)
	 * @throws IOException if it can't rite the files it needs
	 */
	public static void plotFocalPlane(double[] energies, double[][] particles,
									  double[] spectrum,
									  double[] slitPositions, double[] slitLengths,
									  double[] slitWidths) throws IOException {
		CSV.writeColumn(energies, new File("output/focal plane energies.csv"));
		CSV.write(particles, new File("output/focal plane particles.csv"), ',');
		CSV.writeColumn(spectrum, new File("output/focal plane spectrum.csv"));
		CSV.writeColumn(slitLengths, new File("output/focal plane lengths.csv"));
		CSV.writeColumn(slitWidths, new File("output/focal plane widths.csv"));
		CSV.writeColumn(slitPositions, new File("output/focal plane positions.csv"));
		ProcessBuilder plotPB = new ProcessBuilder("python", "src/python/plot_focal_plane.py");
		plotPB.start();
	}
}
