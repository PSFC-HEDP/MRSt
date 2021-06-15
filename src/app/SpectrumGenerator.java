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
package app;

import util.CSV;
import util.NumericalMethods;

import java.io.File;
import java.io.IOException;

/**
 * @author Justin Kunimune
 *
 */
public class SpectrumGenerator {

	private static final NumericalMethods.DiscreteFunction ALPHA_KNOCKON_SPECTRUM = new NumericalMethods.DiscreteFunction(
			new double[] {10.23, 10.5, 11.0, 11.25, 11.5, 12.0, 12.5, 13.0, 13.5, 14.0, 14.5,
					15.0, 15.5, 16.0, 16.5, 17.0, 17.5, 18.0, 18.5, 19.0, 19.5, 19.85},
			new double[] {1.62E-06, 4.87E-06, 3.71E-05, 8.85E-05, 0.00024044, 0.0019635, 0.016034, 0.097, 0.17674, 0.21588, 0.21588,
					0.17674, 0.071859, 0.019584, 0.0056109, 0.00169, 0.00046811, 0.00014583, 4.26E-05, 1.28E-05, 3.68E-06, 1.62E-06}); // [1/MeV]
	private static final NumericalMethods.DiscreteFunction DOWN_SCATTER_SPECTRUM = new NumericalMethods.DiscreteFunction(
			new double[] {11.50, 11.75, 12.00, 12.25, 12.50, 12.75, 13.00, 13.25,
					13.50, 13.75, 14.00, 14.25, 14.50},
			new double[] {0.026877796, 0.029223872, 0.030997082, 0.033544329, 0.035526223, 0.038301112, 0.040480957, 0.043125867,
					0.045434499, 0.048972573, 0.05105225, 0, 0}, true); // [1/MeV/(g/cm^2)]

	/**
	 * generate a time-averaged spectrum based on some parameters that are taken to be constant.
	 * @param Yn the total neutron yield [10^15]
	 * @param Ti the ion temperature [keV]
	 * @param Te the electron temperature [keV]
	 * @param vi the bulk flow rate parallel to the line of sight [μm/ns]
	 * @param ρR the areal density of fuel and shell surrounding the hot spot [g/cm^2]
	 * @param eBins the edges of the energy bins [MeV]
	 * @return the theoretical number of particles in each energy bin, ignoring stochastity.
	 */
	public static double[] generateSpectrum(
			double Yn, double Ti, double Te, double vi, double ρR, double[] eBins) {
		return generateSpectrum(Yn, Ti, Te, vi, ρR, eBins, false, null);
	}

	/**
	 * generate a time-averaged spectrum based on some parameters that are taken to be constant.
	 * @param Yn the total neutron yield [10^15]
	 * @param Ti the ion temperature [keV]
	 * @param Te the electron temperature [keV]
	 * @param vi the bulk flow rate parallel to the line of sight [μm/ns]
	 * @param ρR the areal density of fuel and shell surrounding the hot spot [g/cm^2]
	 * @param eBins the edges of the energy bins [MeV]
	 * @return the theoretical number of particles in each energy bin, ignoring stochastity.
	 */
	public static double[] generateSpectrum(
			double Yn, double Ti, double Te, double vi, double ρR,
			double[] eBins, boolean onlyDS) {
		return generateSpectrum(Yn, Ti, Te, vi, ρR, eBins, onlyDS, null);
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
			double[] timeSlice = generateSpectrum(Yn[j]*dt, Ti[j], Te[j], vi[j], ρR[j], eBins);
			for (int i = 0; i < spectrum.length; i ++)
				spectrum[i][j] = timeSlice[i];
		}
		return spectrum;
	}


	/**
	 * generate a time-averaged spectrum based on some parameters that are taken to be constant.
	 * @param Yn the primary neutron yield [10^15]
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
			double[] eBins, boolean onlyDS, NumericalMethods.DiscreteFunction downScatterCalibration) {
		double ΔEth = 5.30509e-3/(1 + 2.4736e-3*Math.pow(Ti, 1.84))*Math.pow(Ti, 2/3.) + 1.3818e-3*Ti;
		double μ = Math.max(0, 14.029 + ΔEth + .54e-3*vi); // primary peak (see paper) [MeV]
		double σ2 = .4034*μ*Ti/1e3; // primary width [MeV^2]
		double upscat = 1 - Math.exp(-8.6670e-5*Math.pow(Te, 2.5149)); // probability of a neutron being scattered up by an alpha
		double downscat = 1 - Math.exp(-.255184*ρR); // probability of a neutron being scattered down by DT
		double total = Yn/(1 - upscat)/(1 - downscat); // total yield

		double[] I = new double[eBins.length]; // probability distribution at edges [MeV^-1]
		for (int i = 0; i < eBins.length; i ++) {
			if (Ti > 0 && σ2 > 0) {
				if (!onlyDS) {
					I[i] += Yn*1e15/Math.sqrt(2*Math.PI*σ2)*
							Math.exp(-Math.pow((eBins[i] - μ), 2)/(2*σ2));
					I[i] += upscat*total*1e15*ALPHA_KNOCKON_SPECTRUM.evaluate(eBins[i]);
				}
				I[i] += ρR*total*1e15*DOWN_SCATTER_SPECTRUM.evaluate(eBins[i]);
			}
			else
				I[i] = 0;
		}

		double[] counts = new double[eBins.length-1];
		for (int i = 0; i < counts.length; i ++)
			counts[i] = (I[i] + I[i+1])/2.*(eBins[i+1] - eBins[i]);
		return counts;
	}

	/**
	 * convert a time-cumulative spectrum, as read from an input file, into an array of counts.
	 * also, edit energies to make it bin edges. this will remove the last row of energies.
	 * given that the last row is over 29 MeV, I think that's fine.
	 * @param times the time bin boundaries [ns]
	 * @param energies the energy bin midpoints [MeV]
	 * @param spectrum array of the neutron spectrum integrated in time [#/MeV]
	 * @return spectrum array of the neutron spectrum [#]
	 */
	public static double[][]
	interpretSpectrumFile(double[] times, double[] energies, double[][] spectrum) {
		for (int i = energies.length-1; i > 0; i --) // start by fixing the energies
			energies[i] = (energies[i-1] + energies[i])/2;
		energies[0] = 2*energies[0] - energies[1]; // for the last one, assume equal spacing

		double[][] output = new double[energies.length-1][times.length-1];
		for (int i = 0; i < energies.length-1; i ++) { // now for the spectrum
			for (int j = 0; j < times.length-1; j ++)
				output[i][j] = (spectrum[i][j+1] - spectrum[i][j])*(energies[i+1] - energies[i]); // you're just taking this derivative
		}

		return output;
	}

	/**
	 * modify a spectrum artificially in place.
	 * @param times the time bin boundaries [ns]
	 * @param energies the energy bin midpoints [MeV]
	 * @param spectrum array of the neutron spectrum integrated in time [#/MeV]
	 * @param yield the flat yield modifier to apply
	 * @param temp the flat temperature modifier to apply
	 * @param downS the flat down scatter yield modifier to apply
	 * @param flow the flat velocity shift to apply [km/s]
	 */
	public static void modifySpectrum(double[] times, double[] energies, double[][] spectrum,
	                                  double yield, double temp, double downS, double flow) {
		for (int i = 0; i < energies.length - 1; i ++) { // scale the whole thing up or down to change yield (and account for the broadening)
			for (int j = 0; j < times.length - 1; j ++) {
				if (energies[i] >= 13.3)
					spectrum[i][j] = yield/Math.sqrt(temp)*spectrum[i][j];
				else
					spectrum[i][j] = yield*downS/Math.sqrt(temp)*spectrum[i][j];
			}
		}
		for (int j = 0; j < times.length - 1; j ++) { // scale it in energy space to change temperature
			double[] slice = new double[energies.length - 1];
			for (int i = 0; i < energies.length - 1; i ++)
				slice[i] = spectrum[i][j];
			int argmax = NumericalMethods.argmax(slice);
			double ePeak = NumericalMethods.mean(energies, slice);
			double iPeak = argmax + (ePeak - energies[argmax])/(energies[1] - energies[0]); // find the peak index (assume equally spaced energy bins)
			iPeak = NumericalMethods.coerce(argmax-1, argmax+1, iPeak);
			for (int i = 0; i < energies.length - 1; i ++) {
				double iP = iPeak + (i - iPeak)/Math.sqrt(temp);
				if (iP >= 0 && iP < energies.length-2)
					spectrum[i][j] = (1-iP%1)*slice[(int)iP] + (iP%1)*slice[(int)iP+1];
				else
					spectrum[i][j] = 0;
			}
		}
		for (int i = 0; i < energies.length; i ++)
			energies[i] = energies[i] - .54e-3*flow;
	}


	/**
	 * Generate the spectra files needed for analysis from the corresponding
	 * trajectory CSVs.
	 */
	public static void main(String[] args) throws NumberFormatException, IOException {
		for (String filename : new String[] {"failed", "marginal", "high", "og", "og with falling temp"}) {
			double[][] thing;
			double[] eBins;
			thing = CSV.read(new File("input/trajectories "+filename+".csv"), ',', 1);
			eBins = CSV.readColumn(new File("input/energy.txt"));
			
			double[] time = new double[thing.length];
			double[] ρR = new double[thing.length];
			double[] Yn = new double[thing.length];
			double[] Ti = new double[thing.length];
			double[] zero = new double[thing.length];
			for (int i = 0; i < thing.length; i ++) {
				time[i] = thing[i][0];
				Yn[i] = thing[i][1]*.1*1e6/1e-6/(14e6*1.6e-19)/1e15*1e-9; // convert from 0.1MJ/μs to 1e15n/ns
				Ti[i] = thing[i][4];
				ρR[i] = thing[i][3];
				zero[i] = 0;
			}
			double[] tBins = new double[time.length + 1];
			tBins[0] = (3*time[0] - time[1])/2.;
			for (int i = 1; i < time.length; i ++)
				tBins[i] = (time[i-1] + time[i])/2.;
			tBins[time.length] = (3*time[time.length-1] - time[time.length-2])/2.;
			double[][] spectrum = generateSpectrum(Yn, Ti, zero, zero, ρR, eBins, tBins);
			
			CSV.writeColumn(tBins, new File("input/time "+filename+".txt"));
			CSV.write(spectrum, new File("input/spectrum "+filename+".txt"), '\t');
		}
		
		System.out.println("done");
	}

}
