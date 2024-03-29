/**
 * MIT License
 *
 * Copyright (c) 2020 Justin Kunimune
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
package physics;

import util.CSV;
import util.Math2;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

/**
 * @author Justin Kunimune
 *
 */
public class SpectrumGenerator {

	private static final Math2.DiscreteFunction ALPHA_KNOCKON_SPECTRUM = new Math2.DiscreteFunction(
			new double[] {10.23, 10.5, 11.0, 11.25, 11.5, 12.0, 12.5, 13.0, 13.5, 14.0, 14.5,
					15.0, 15.5, 16.0, 16.5, 17.0, 17.5, 18.0, 18.5, 19.0, 19.5, 19.85},
			new double[] {1.62E-06, 4.87E-06, 3.71E-05, 8.85E-05, 0.00024044, 0.0019635, 0.016034, 0.097, 0.17674, 0.21588, 0.21588,
					0.17674, 0.071859, 0.019584, 0.0056109, 0.00169, 0.00046811, 0.00014583, 4.26E-05, 1.28E-05, 3.68E-06, 1.62E-06},
			false, true); // [1/MeV]
	private static final Math2.DiscreteFunction DOWN_SCATTER_SPECTRUM = new Math2.DiscreteFunction(
			new double[] {11.50, 11.75, 12.00, 12.25, 12.50, 12.75, 13.00, 13.25,
					13.50, 13.75, 14.00, 14.25, 14.50},
			new double[] {0.026877796, 0.029223872, 0.030997082, 0.033544329, 0.035526223, 0.038301112, 0.040480957, 0.043125867,
					0.045434499, 0.048972573, 0.05105225, 0, 0},
			true, false); // [1/MeV/(g/cm^2)]


	/**
	 * calculate the energy at which the downscatter spectrum overtakes the primary
	 * @param Ti the ion temperature (keV)
	 * @param ρR the areal density (g/cm^2)
	 * @return the point at which the two components are equal, accurate to 10 keV
	 */
	public static double primaryCutoff(double Ti, double ρR) {
		double[] energy = new double[201];
		for (int i = 0; i < energy.length; i ++)
			energy[i] = 12. + 2.*i/(energy.length - 1);
		double[] total = generateSpectrum(1, Ti, 0, 0, ρR, energy, false);
		double[] secondary = generateSpectrum(1, Ti, 0, 0, ρR, energy, true);
		for (int i = 0; i < energy.length - 1; i ++)
			if (total[i] > 2*secondary[i])
				return energy[i];
		throw new RuntimeException("this is a nonphysically hi ρR");
	}

	/**
	 * generate a time-averaged spectrum based on some parameters that are taken to be constant.
	 * @param Yn the primary neutron yield []
	 * @param Ti the ion temperature [keV]
	 * @param Te the electron temperature [keV]
	 * @param vi the bulk flow rate parallel to the line of sight [μm/ns]
	 * @param ρR the areal density of fuel and shell surrounding the hot spot [g/cm^2]
	 * @param eBins the edges of the energy bins [MeV]
	 * @return the theoretical number of particles in each energy bin, ignoring stochastity.
	 */
	public static double[] generateSpectrum(
			double Yn, double Ti, double Te, double vi, double ρR, double[] eBins) {
		return generateSpectrum(Yn, Ti, Te, vi, ρR, eBins, false);
	}


	/**
	 * generate a time-resolved spectrum based on some parameters that vary with time.
	 * @param Yn the neutron yield rate [10^15/ns]
	 * @param Ti the ion temperature [keV]
	 * @param Te the electron temperature [keV]
	 * @param vi the bulk flow rate parallel to the line of sight [μm/ns]
	 * @param ρR the areal density of fuel and shell surrounding the hot spot [g/cm^2]
	 * @param tBins the edges of the time bins [ns]
	 * @param eBins the edges of the energy bins [MeV]
	 * @return the theoretical number of particles in each energy bin, ignoring stochastity.
	 */
	public static double[][] generateSpectrum(
			double[] Yn, double[] Ti, double[] Te, double[] vi, double[] ρR,
			double[] eBins, double[] tBins) {
		double[][] spectrum = new double[eBins.length-1][tBins.length-1];
		for (int j = 0; j < spectrum[0].length; j ++) {
			double dt = (tBins[j+1] - tBins[j]); // bin width [ns]
			double[] timeSlice = generateSpectrum(Yn[j]*1e15*dt, Ti[j], Te[j], vi[j], ρR[j], eBins);
			for (int i = 0; i < spectrum.length; i ++)
				spectrum[i][j] = timeSlice[i];
		}
		return spectrum;
	}


	/**
	 * generate a time-averaged spectrum based on some parameters that are taken to be constant.
	 * @param Yn the primary neutron yield []
	 * @param Ti the ion temperature [keV]
	 * @param Te the electron temperature [keV]
	 * @param vi the bulk flow rate parallel to the line of sight [μm/ns]
	 * @param ρR the areal density of fuel and shell surrounding the hot spot [g/cm^2]
	 * @param eBins the edges of the energy bins [MeV]
	 * @param onlyDS only return the DS spectrum; remove all primaries
	 * @return the theoretical number of particles in each energy bin, ignoring stochastity.
	 */
	public static double[] generateSpectrum(
			double Yn, double Ti, double Te, double vi, double ρR,
			double[] eBins, boolean onlyDS) {
		double ΔEth = 5.30509e-3/(1 + 2.4736e-3*Math.pow(Ti, 1.84))*Math.pow(Ti, 2/3.) + 1.3818e-3*Ti;
		double μ = Math.max(0, 14.029 + ΔEth + .54e-3*vi); // primary peak (see paper) [MeV]
		double σ2 = Math.abs(.4034*μ*Ti/1e3); // primary width [MeV^2]
		double downscat = 1 - Math.exp(-.2107*ρR); // probability of a neutron being scattered down by DT
		double upscat = (1 - Math.exp(-8.6670e-5*Math.pow(Te, 2.5149)))*(1 - downscat); // probability of a neutron being scattered up by an alpha
		double primary = 1 - upscat - downscat; // probability of a neutron escaping unscatterd
		double total = Yn/Math.max(primary, 0.1); // total yield (here's a weerd edge case: if the ρR is huge, don’t correct this by more than ×10)

		double[] I = new double[eBins.length]; // calculate spectrum density at the edges
		double[] IPrimary = new double[eBins.length];
		for (int i = 0; i < eBins.length; i ++) {
			double primaryComponent = (σ2 > 0) ?
					total*primary/Math.sqrt(2*Math.PI*σ2)*
							Math.exp(-Math.pow((eBins[i] - μ), 2)/(2*σ2)) : 0;
			double downscatComponent = total*ρR*DOWN_SCATTER_SPECTRUM.evaluate(eBins[i]);
			double upscatComponent = total*upscat*ALPHA_KNOCKON_SPECTRUM.evaluate(eBins[i]);
			if (onlyDS)
				I[i] = downscatComponent;
			else
				I[i] = downscatComponent + primaryComponent + upscatComponent;
			IPrimary[i] = primaryComponent;
		}

		double primaryTotal = 0;
		double[] counts = new double[eBins.length - 1];
		for (int i = 0; i < eBins.length - 1; i ++) {
			counts[i] = (I[i] + I[i+1])/2*(eBins[i+1] - eBins[i]);
			primaryTotal += (IPrimary[i] + IPrimary[i+1])/2*(eBins[i+1] - eBins[i]);
		}
		int peakBin = Math2.bin(μ, eBins);
		if (primaryTotal < total*primary && peakBin >= 0)
			counts[peakBin] += total*primary - primaryTotal; // make up for any particles lost to curvature

		for (double value: counts)
			if (Double.isNaN(value))
				throw new IllegalArgumentException(String.format("passing Yn=%.4g, Ti=%.4g, Te=%.4g, vi=%.4g, rhoR=%.4g yields a nan.", Yn, Ti, Te, vi, ρR));
		return counts;
	}

	/**
	 * check the spectrum to see if it is time-cumulative. if it is, convert it
	 * into an array of counts. also, check the energies to see if they are bin
	 * centers, and if so, edit energies to make it bin edges. this will remove
	 * the last row of energies. given that the last row is over 29 MeV, I think
	 * that's fine.
	 * @param times the time bin boundaries [ns]
	 * @param energies the energy bin midpoints [MeV]
	 * @param spectrum array of the neutron spectrum integrated in time (or not) [#/MeV]
	 * @return spectrum array of the neutron spectrum [#]
	 */
	public static double[][] interpretSpectrumFile(
			double[] times, double[] energies, double[][] spectrum) {
		if (spectrum.length != energies.length - 1) {
			if (spectrum.length == energies.length) {
				System.err.println("WARNING: these energies were bin centers. why?");
				for (int i = energies.length - 1; i > 0; i--) // start by fixing the energies
					energies[i] = (energies[i - 1] + energies[i])/2;
				energies[0] = 2*energies[0] - energies[1]; // for the last one, assume equal spacing
			}
			else {
				throw new IllegalArgumentException("the spectrum shape "+spectrum.length+"x"+spectrum[0].length+" is incompatible with "+energies.length+" energies");
			}
		}
		if (spectrum[0].length != times.length - 1) {
			if (spectrum[0].length == times.length) {
				System.err.println("WARNING: this spectrum was cumulative. why??");
				double[][] output = new double[energies.length - 1][times.length - 1];
				for (int i = 0; i < energies.length - 1; i++) { // now for the spectrum
					for (int j = 0; j < times.length - 1; j++)
						output[i][j] = (spectrum[i][j + 1] - spectrum[i][j])*(energies[i + 1] - energies[i]);
				}
				return output;
			}
			else {
				throw new IllegalArgumentException("the spectrum shape "+spectrum.length+"x"+spectrum[0].length+" is incompatible with "+times.length+" times");
			}
		}
		else {
			return spectrum;
		}
	}

	/**
	 * modify a spectrum artificially in place.
	 * @param spectrum array of the neutron spectrum integrated in time [#/MeV]
	 * @param yieldAmplification the amount to scale it
	 */
	public static double[][] modifySpectrum(double[][] spectrum,
	                                        double yieldAmplification) {
		double[][] output = new double[spectrum.length][spectrum[0].length];
		for (int i = 0; i < spectrum.length; i ++) // scale the whole thing up or down to change yield
			for (int j = 0; j < spectrum[i].length; j ++)
				output[i][j] = yieldAmplification*spectrum[i][j];
		return output;
	}


	/**
	 * Generate the spectra files needed for analysis from the corresponding
	 * trajectory CSVs.
	 */
	public static void main(String[] args) throws NumberFormatException, IOException {
		//		double[] eBins = CSV.readColumn(new File("input/energy.txt"));
		//		double[] eBins = new double[(int) ((16 - 12)/.05 + 1)];
		//		for (int i = 0; i < eBins.length; i ++)
		//			eBins[i] = (12 + i*(16. - 12.)/(eBins.length-1));
		//		double[] primary = generateSpectrum(5.2e21, 3, 4, 0, 0.0, eBins);
		//		double[] full = generateSpectrum(1, 3, 4, 0, 3.0, eBins);
		//		for (int i = 0; i < eBins.length - 1; i ++) {
		//			System.out.printf("[%f, %f, %f],\n", (eBins[i] + eBins[i+1])/2, primary[i], full[i]);
		//		}
		for (Path file : (Iterable<Path>) Files.walk(Paths.get("input/"))::iterator) {
			if (file.getFileName().toString().contains("trajectories")) {
				System.out.println(file.getFileName());

				String[] header = CSV.readHeader(file.toFile(), ',');
				double[][] thing;
				try {
					thing = CSV.read(file.toFile(), ',', 1);
				} catch (NumberFormatException e) {
					thing = CSV.read(file.toFile(), ' ', 1);
				}
				double[] eBins = CSV.readColumn(new File("input/energy.txt"));

				int timeIndex = -1, burnIndex = -1, tionIndex = -1, totalDensityIndex = -1;
				for (int i = 0; i < header.length; i ++) {
					String head = header[i].toLowerCase();
					if (head.contains("time"))
						timeIndex = i;
					else if (head.contains("burn") || head.contains("eprod") || head.contains("neut"))
						burnIndex = i;
					else if (head.contains("tion"))
						tionIndex = i;
					else if ((head.contains("total") && head.contains("rhor")) || head.contains("dsr"))
						totalDensityIndex = i;
				}
				double timeUnit;
				if (header[timeIndex].contains("us"))
					timeUnit = 1e-6;
				else if (header[timeIndex].contains("ns"))
					timeUnit = 1e-9;
				else
					timeUnit = 1e-15;
				double burnUnit;
				if (header[burnIndex].contains("us^-1") || header[burnIndex].contains("per us"))
					burnUnit = 1e+6;
				else if (header[burnIndex].contains("ns^-1") || header[burnIndex].contains("per ns"))
					burnUnit = 1e+9;
				else if (header[burnIndex].contains("ps^-1") || header[burnIndex].contains("per ps"))
					burnUnit = 1e+12;
				else if (header[burnIndex].contains("kJ/ps") || header[timeIndex].contains("MJ/ns"))
					burnUnit = 1e+15/(14e6*1.6e-19);
				else // 0.1MJ/μs
					burnUnit = .1*1e+12/(14e6*1.6e-19);
				double ρRUnit;
				if (header[totalDensityIndex].toLowerCase().contains("dsr"))
					ρRUnit = 20.4;
				else if (header[totalDensityIndex].contains("mg"))
					ρRUnit = 1e-3;
				else
					ρRUnit = 1e-0;

				double[] time = new double[thing.length];
				double[] ρR = new double[thing.length];
				double[] Yn = new double[thing.length];
				double[] Ti = new double[thing.length];
				double[] zero = new double[thing.length];
				for (int i = 0; i < thing.length; i++) {
					time[i] = thing[i][timeIndex]*timeUnit/1e-9; // convert to ns
					Yn[i] = thing[i][burnIndex]*burnUnit/(1e15/1e-9); // convert to (1e15/ns)
					Ti[i] = thing[i][tionIndex];
					ρR[i] = thing[i][totalDensityIndex]*ρRUnit; // convert to g/cm^2
					System.out.println(ρR[i]);
					zero[i] = 0;
				}
				double[] tBins = new double[time.length + 1];
				tBins[0] = (3*time[0] - time[1])/2.;
				for (int i = 1; i < time.length; i++)
					tBins[i] = (time[i - 1] + time[i])/2.;
				tBins[time.length] = (3*time[time.length - 1] - time[time.length - 2])/2.;

				double[][] spectrum = generateSpectrum(Yn, Ti, zero, zero, ρR, eBins, tBins);

				CSV.writeColumn(tBins, new File(file.toString()
						.replace("trajectories", "time").replace(".csv", ".txt")));
				CSV.write(spectrum, new File(file.toString()
						.replace("trajectories", "spectrum").replace(".csv", ".txt")), '\t');
			}
		}

		System.out.println("done");
	}

}
