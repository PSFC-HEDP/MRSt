/**
 * MIT License
 * 
 * Copyright (c) 2018 Justin Kunimune
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
package main;

import java.io.IOException;
import java.util.Arrays;
import java.util.Locale;
import java.util.logging.Logger;

import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.MaxIter;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.MultivariateOptimizer;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.MultiDirectionalSimplex;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.PowellOptimizer;

import main.NumericalMethods.DiscreteFunction;

/**
 * the class where all the math is.
 * 
 * @author Justin Kunimune
 */
public class MRSt {
	
	private static final int x = 0, y = 1, z = 2;
	
	private static final double SPEED_OF_LIGHT = 2.99792458e8;
	
	private static final int STOPPING_DISTANCE_RESOLUTION = 64;
	private static final double MIN_E = 12, MAX_E = 16; // histogram bounds [MeV]
	private static final double MIN_T = 16.0, MAX_T = 16.5; // histogram bounds [ns]
	private static final double E_RESOLUTION = .1, T_RESOLUTION = 13e-3; // resolutions [MeV], [ns]
//	private static final double E_RESOLUTION = .3, T_RESOLUTION = 40e-3;
	private static final int MIN_STATISTICS = 100; // the minimum number of deuterons to define a spectrum at a time
	private static final int TRANSFER_MATRIX_TRIES = 10000; // the number of points to sample in each column of the transfer matrix
	
	private static final double[] PARAM_SCALES = { Double.NaN, 4, 100, 1 }; // I'm 68% sure the magnitudes of Y, v, T, and ρR won't exceed these
	
	private final double foilDistance; // z coordinate of midplane of foil [m]
	private final double foilThickness; // thickness of foil [m]
	private final double apertureDistance; // distance from TCC to aperture [m]
	private final double apertureWidth; // horizontal dimension of aperture [m]
	private final double apertureHeight; // vertical dimension of aperture [m]
	private final double cosyKmin, cosyKmax; // bounds on deuterons that will be detected by the CSI
	private final double cosyK0, cosyV0, cosyT0, cosyT1; // reference coordinates for beam used by COSY [J], [m/s], [s], [s?]
	private final double focalPlaneAngle; // angle of focal plane [rad]
	private final double[][] cosyCoefficients; // table of coefficients generated by COSY
	private final int[][] cosyExponents; // table of term powers generated by COSY
	private final Particle ion; // either D or P, but probably D

	private final double energyFactor; // conversion factor between neutron and ion energies []
	private final double probHitsFoil; // probability that the neutron goes through the foil
	
	private final DiscreteFunction distanceVsEnergy; // stopping distance info
	private final DiscreteFunction energyVsDistance; // inverse stopping distance info
	private final DiscreteFunction energyVsPosition; // map between location on detector and energy going into lens
	
	private final double[] energyBins; // endpoints of E bins for inferred spectrum [MeV]
	private final double[] timeBins; // endpoints of time bins for inferred spectrum [ns]
	private final double[][] transferMatrix; // the full nmx2nm transfer matrix plus smoothing rows
	private double[][] deuteronSpectrum; // time-corrected deuteron counts
	private double[][] fitNeutronSpectrum; // backward-fit neutron counts
	private double[][] fitDeuteronSpectrum; // backward-fit deuteron counts (this should be similar to deuteronSpectrum)
	
	private final double[] timeAxis; // 1D vectors for higher level measurements [ns]
	private double[] ionTemperature; // [keV]
	private double[] flowVelocity; // [km/s]
	private double[] arealDensity; // [g/cm^2]
	private double[] neutronYield; // [10^15/ns]
//	private double[] neutronMean; // [keV]
//	private double[] neutronWidth; // [keV^2]
//	private double[] neutronSkew; // [keV^3]
//	private double[] neutronKurtosis; // [keV^4]
	
	private final Logger logger; // for logging
	
	/**
	 * perform some preliminary calculations for the provided configuration.
	 * @param stoppingPower a double[][] containing two columns and n rows. the zeroth column is
	 * the reference values of E in [keV] and the last column is the corresponding values of
	 * dE/dx in [keV/μm].
	 * @param focalTilt angle of the focal plane (0 means untilted) [deg]
	 * @param referenceEnergy minimum ion energy allowed in COSY calculation [eV]
	 * @param referenceEnergy maximum ion energy allowed in COSY calculation [eV]
	 * @param referenceEnergy expected ion energy used in COSY calculation [eV]
	 */
	public MRSt(
			Particle ion, double foilDistance, double foilRadius, double foilThickness,
			double[][] stoppingPowerData,
			double apertureDistance, double apertureWidth, double apertureHeight,
			double minimumEnergy, double maximumEnergy, double referenceEnergy,
			double[][] cosyCoefficients, int[][] cosyExponents,
			double focalTilt, int numBins, Logger logger) {
		this.foilDistance = foilDistance;
		this.foilThickness = foilThickness;
		this.apertureDistance = apertureDistance;
		this.apertureWidth = apertureWidth;
		this.apertureHeight = apertureHeight;
		this.cosyKmin = minimumEnergy*(-Particle.E.charge);
		this.cosyKmax = maximumEnergy*(-Particle.E.charge);
		this.cosyK0 = referenceEnergy*(-Particle.E.charge); // save this in a more useful unit
		this.cosyV0 = Math.sqrt(2*cosyK0/ion.mass); // and get the corresponding speed
		double γ = Math.pow(1 - Math.pow(cosyV0/SPEED_OF_LIGHT, 2), -1/2.);
		double L = (1 + γ)/γ*2*cosyCoefficients[5][4]; // here's a fun shortcut to estimating the length of the lens: first order analysis
		this.cosyT0 = L/cosyV0; // and corresponding time
		this.cosyT1 = -(1+γ)/γ/cosyV0; // and why is time measured in units of distance?
		this.focalPlaneAngle = Math.toRadians(focalTilt);
		this.cosyCoefficients = cosyCoefficients;
		this.cosyExponents = cosyExponents;
		this.ion = ion;
		
		double A = ion.mass/Particle.N.mass;
		this.energyFactor = 4*A/Math.pow(A + 1, 2);
		
		double foilMaxAngle = Math.atan(foilRadius/foilDistance);
		this.probHitsFoil = (1 - Math.cos(foilMaxAngle))/2;
		
		double[] dxdE = new double[stoppingPowerData.length]; // integrate the stopping power to get stopping distance
		double[] E = new double[stoppingPowerData.length];
		for (int i = 0; i < stoppingPowerData.length; i ++) {
			dxdE[i] = 1/(stoppingPowerData[i][1]*1e9*(-Particle.E.charge)); // converting from [GeV/m] to [m/J]
			E[i] = stoppingPowerData[i][0]*1e3*(-Particle.E.charge); // and from [keV] to [J]
		}
		DiscreteFunction distanceVsEnergyRaw = new DiscreteFunction(E, dxdE).antiderivative();
		this.distanceVsEnergy = distanceVsEnergyRaw.indexed(STOPPING_DISTANCE_RESOLUTION); // m(J)
		this.energyVsDistance = distanceVsEnergyRaw.inv().indexed(STOPPING_DISTANCE_RESOLUTION); // J(m)
		
		this.energyBins = new double[(int) ((MAX_E - MIN_E)/E_RESOLUTION + 1)];
		for (int i = 0; i < energyBins.length; i ++)
			this.energyBins[i] = (MIN_E + i*(MAX_E - MIN_E)/(energyBins.length-1));
		this.timeBins = new double[(int) ((MAX_T - MIN_T)/T_RESOLUTION + 1)];
		for (int i = 0; i < timeBins.length; i ++)
			this.timeBins[i] = MIN_T + i*(MAX_T - MIN_T)/(timeBins.length-1);
		
		this.timeAxis = new double[timeBins.length-1];
		for (int i = 0; i < timeBins.length-1; i ++)
			this.timeAxis[i] = (this.timeBins[i] + this.timeBins[i+1])/2;
		
		double[] calibEnergies = new double[2*energyBins.length];
		double[] detectorPosition = new double[calibEnergies.length];
		for (int i = 0; i < calibEnergies.length; i ++) {
			calibEnergies[i] = (MIN_E + (MAX_E - MIN_E)*i/(calibEnergies.length-1))*energyFactor*1e6*(-Particle.E.charge);
			double[] v = {0, 0, Math.sqrt(2*calibEnergies[i]/ion.mass)};
			double[] r = computeFocusedPosition(new double[] {0,0,0}, v, 0);
			detectorPosition[i] = r[x]/Math.cos(focalPlaneAngle);
		}
		this.energyVsPosition = new DiscreteFunction(calibEnergies, detectorPosition).inv()
				.indexed(energyBins.length); // J(m)
		
		this.logger = logger;
		
		this.transferMatrix = evaluateTransferMatrix();
	}
	
	/**
	 * compute the probability that a given neutron released at this energy will spawn an ion
	 * and knock that ion through the aperture.
	 * @param energy energy of the released particles [eV]
	 * @return the fraction of particles that are worth simulating
	 */
	public double efficiency(double energy) {
		double n = 0.08e2; // I'm not sure what units this has or whence it came
		double dσdΩ = 4.3228e3/Math.sqrt(energy) - 0.6523; // same with these ones
		double dΩ = apertureWidth*apertureHeight / Math.pow(apertureDistance - foilDistance, 2);
		double l = foilThickness; // assume the foil is thin so we don't have to worry about multiple collisions
		return probHitsFoil * n*dσdΩ*dΩ*l;
	}
	
	/**
	 * determine the response of the MRSt to neutrons at particular energies. also tack on some
	 * finite difference laplacian roughness measures.
	 * @return the matrix that multiplies from a real neutron spectrum to a measured deuteron
	 *   spectrum concatenate with roughness, both flattened out.
	 */
	public double[][] evaluateTransferMatrix() {
		if (logger != null)  logger.info("beginning Monte Carlo computation.");
		int simulationCount = 0;
		long startTime = System.currentTimeMillis();
		
		int n = (energyBins.length-1)*(timeBins.length-1);
		double[][] matrix = new double[n][n];
		int j0 = timeBins.length/2; // time symmetry means we only need to evaluate at one time
		for (int i = 0; i < energyBins.length-1; i ++) { // sweep through all energies
			double energy0 = energyBins[i]*1e6, energy1 = energyBins[i+1]*1e6; // [eV]
			double time0 = timeBins[j0]*1e-9, time1 = timeBins[j0+1]*1e-9; // [s]
			double weight = efficiency((energy0 + energy1)/2)/TRANSFER_MATRIX_TRIES;
			for (int k = 0; k < TRANSFER_MATRIX_TRIES; k ++) {
				double energy = energy0 + Math.random()*(energy1 - energy0); // randomly choose values from the bin
				double time = time0 + Math.random()*(time1 - time0); // [s]
				
				double[] xt = this.simulate(energy, time); // do the simulation!
				simulationCount ++;
				
				if (Double.isNaN(xt[0]))	continue; // sometimes, they won't hit the CSI. That's fine.
				double[] et = this.backCalculate(xt[0], xt[1]); // otherwise do the stretch/compress time correction
				
				double e = et[0]/(-Particle.E.charge)/1e6, t = et[1]/1e-9; // then convert to the same units as the bins
//				e = energy/1e6; t = time/1e-9;
				int eBin = (int)((e - MIN_E)/(MAX_E - MIN_E)*(energyBins.length-1));
				int tBin = (int)((t - MIN_T)/(MAX_T - MIN_T)*(timeBins.length-1));
				if (eBin >= 0 && eBin < energyBins.length-1 && tBin >= 0 && tBin < timeBins.length-1) // if it falls in detectable bounds
					matrix[(timeBins.length-1)*eBin + tBin][(timeBins.length-1)*i + j0] += weight; // add it to the row
			}
		}
		
		long endTime = System.currentTimeMillis();
		if (logger != null)
			logger.info(String.format(Locale.US, "completed %d simulations in %.2f minutes.",
					simulationCount, (endTime - startTime)/60000.));
		
		for (int i = 0; i < n; i ++) { // now iterate through all of the rows
			for (int j = 0; j < n; j ++) { // and all of the columns, of which we still don't have many
				int jRef = j/(timeBins.length-1)*(timeBins.length-1) + j0; // look for the nearest column that is filled
				int iRef = i - j + jRef; // [i,j] is equivalent to an [iRef,jRef] by time symmetry
				if (iRef >= 0 && i/(timeBins.length-1) == iRef/(timeBins.length-1)) // this may have jumped over to the next energy,
					matrix[i][j] = matrix[iRef][jRef];
				else
					matrix[i][j] = 0; // but if so, it's almost certainly 0
			}
		}
		
		return matrix;
	}

	/**
	 * compute the total response to an implosion with the given neutron spectrum and save and
	 * analyze it.
	 * @param energies the energies that describe the rows of counts [MeV]
	 * @param times the times that describe the columns of counts [ns]
	 * @param spectrum the time- and energy- resolved neutron spectrum in number of neutrons. each
	 * row corresponds to one element of energies, and each column one element of times. [#/MeV/ns]
	 * @return {computation time, BT, peak-ρR, peak-ρR-ramp, Ti(BT), ρR(BT), vi(BT),
	 *   Ti-ramp(BT), ρR-ramp(BT), vi-ramp(BT), max ρR, yield, μ1, μ2, μ3}
	 */
	public double[] respond(double[] energies, double[] times, double[][] spectrum) {
		this.deuteronSpectrum = this.response(energies, times, spectrum, true);
		return analyze(deuteronSpectrum);
	}
	
	/**
	 * compute the total response to an implosion with the given neutron spectrum using the
	 * precomputed transfer matrix. account for electrostatic time correction, but not for any
	 * analysis.
	 * @param energies the edges of the energy bins
	 * @param times the edges of the time bins
	 * @param neutronSpectrum the counts in each bin
	 * @param stochastic whether to add noise to mimic real data
	 * @return the counts in the measured spectrum bins
	 */
	public double[][] response(double[] inEBins, double[] inTBins, double[][] inSpectrum,
			boolean stochastic) {
		if (inSpectrum.length != inEBins.length-1 || inSpectrum[0].length != inTBins.length-1)
			throw new IllegalArgumentException("These dimensions don't make any sense.");
		
		double[][] sizedSpectrum = NumericalMethods.downsample(
				inTBins, inEBins, inSpectrum, this.timeBins, this.energyBins);
		
		double[] u = new double[(energyBins.length-1)*(timeBins.length-1)];
		for (int i = 0; i < energyBins.length-1; i ++)
			for (int j = 0; j < timeBins.length-1; j ++)
				u[(timeBins.length-1)*i + j] = sizedSpectrum[i][j]; // now flatten the spectrum
		
		double[] v = NumericalMethods.matmul(transferMatrix, u); // then do the multiplication
		
		double[][] outSpectrum = new double[energyBins.length-1][timeBins.length-1];
		for (int i = 0; i < energyBins.length-1; i ++)
			for (int j = 0; j < timeBins.length-1; j ++)
				outSpectrum[i][j] = v[(timeBins.length-1)*i + j]; // finally, unflatten. inflate. rounden.
		
		if (stochastic) { // to simulate stochasticity
			for (int i = 0; i < energyBins.length-1; i ++)
				for (int j = 0; j < timeBins.length-1; j ++)
					outSpectrum[i][j] = NumericalMethods.poisson(outSpectrum[i][j]); // just jitter every cell
		}
		
		return outSpectrum;
	}
	
	/**
	 * take the time-corrected spectrum and use it to compute and store time-resolved values
	 * for ion temperature, areal density, and yield.
	 * @param eBins the edges of the energy bins
	 * @param tBins the edges of the time bins
	 * @param spectrum the time-corrected spectrum we want to understand
	 * @return {computation time, BT, peak-ρR, peak-ρR-ramp, Ti(BT), ρR(BT), vi(BT),
	 *   Ti-ramp(BT), ρR-ramp(BT), vi-ramp(BT), max ρR, yield, μ1, μ2, μ3}
	 */
	private double[] analyze(double[][] spectrum) {
		if (spectrum.length != energyBins.length-1 || spectrum[0].length != timeBins.length-1)
			throw new IllegalArgumentException("These dimensions are wrong.");
		
		double[] initialGuess = new double[4*timeAxis.length]; // first we need a half decent guess
		double totalYield = 0; // also the rough mean estimate of neutron spectrum magnitude for entropy purposes
		final double[] noiseVar = new double[timeAxis.length]; // assume some unknown distribution of this variance, so don't freak out when you see n-knockons
		for (int i = 0; i < timeAxis.length; i ++) {
			double[] timeSlice = new double[energyBins.length-1]; // for that we'll want to do some rudimentary analysis
			for (int j = 0; j < timeSlice.length; j ++)
				timeSlice[j] = spectrum[j][i];
			double yieldEstimate = NumericalMethods.definiteIntegral(energyBins, timeSlice)/
					this.efficiency(14.1e6); // measure the yield in the most basic way possible
			totalYield += yieldEstimate; // [#]
			noiseVar[i] = Math.pow(NumericalMethods.max(timeSlice)/1e3, 2);
			initialGuess[i] = yieldEstimate/1e15/(timeBins[i+1] - timeBins[i]); // [10^15/ns]
			for (int k = 1; k < 4; k ++)
				initialGuess[i+k*timeAxis.length] = PARAM_SCALES[k]; // the initial guess for Ti, vi, and ρR can be uniform
		}
		final double spectrumScale = totalYield/(energyBins.length-1)/(timeBins.length-1);
		PARAM_SCALES[0] = totalYield/1e15/(timeBins[1] - timeBins[0])/(timeBins.length-1);
		for (int i = 0; i < timeAxis.length; i ++) {
			initialGuess[i+1*timeAxis.length] = PARAM_SCALES[1]*(1 + initialGuess[i]/PARAM_SCALES[0])/2;
			initialGuess[i+2*timeAxis.length] = 0;
			initialGuess[i+3*timeAxis.length] = PARAM_SCALES[3]*(1 + initialGuess[i]/PARAM_SCALES[0])/2;
		}
		
		double[] lowerBounds = new double[4*timeAxis.length];
		double[] dimensionScale = new double[4*timeAxis.length];
		double[] upperBounds = new double[4*timeAxis.length];
		for (int i = 0; i < timeAxis.length; i ++) {
			lowerBounds[i+0*timeAxis.length] = 0;
			lowerBounds[i+1*timeAxis.length] = 0;
			lowerBounds[i+2*timeAxis.length] = Double.NEGATIVE_INFINITY;
			lowerBounds[i+3*timeAxis.length] = 0;
			for (int k = 0; k < 4; k ++) {
				dimensionScale[i+k*timeAxis.length] = PARAM_SCALES[k]/6;
				upperBounds[i+k*timeAxis.length] = Double.POSITIVE_INFINITY;
			}
		}
		
		if (logger != null)  logger.info("beginning fit process.");
		System.out.println(Arrays.toString(initialGuess)+",");
		long startTime = System.currentTimeMillis();
		MultivariateOptimizer optimizer = new PowellOptimizer(1e-14, 1);
		double[] opt = optimizer.optimize(
				new InitialGuess(initialGuess),
				new MultiDirectionalSimplex(dimensionScale),
				new MaxIter(100000),
				new MaxEval(1000000),
				GoalType.MINIMIZE,
//				new SimpleBounds(lowerBounds, upperBounds),
				new ObjectiveFunction((double[] guess) -> {
//		double[] opt = Optimization.minimizeLBFGSB(
//				(double[] guess) -> {
					if (Math.random() < 3e-4)
						System.out.println(Arrays.toString(guess)+",");
					
					double[][] params = new double[4][timeAxis.length];
					for (int k = 0; k < 4; k ++) // first unpack the state vector
						for (int i = 0; i < timeAxis.length; i ++)
							params[k][i] = guess[i+k*timeAxis.length];
					
					for (int i = 0; i < timeAxis.length; i ++) { // check for illegal (prior = 0) values
						if (params[0][i] < 0)  return Double.POSITIVE_INFINITY;
						if (params[1][i] < 0)  return Double.POSITIVE_INFINITY;
						if (params[3][i] < 0)  return Double.POSITIVE_INFINITY;
					}
					
//					System.out.println(Arrays.toString(guess)+",");
					double[][] teoSpectrum = generateSpectrum(
							params[0], params[1], params[2], params[3], energyBins, timeBins); // generate the neutron spectrum based on those
					double[][] fitSpectrum = this.response(
							energyBins, timeBins, teoSpectrum, false); // blur it according to the transfer matrix
					
					double err = 0; // negative Bayes factor (in nepers)
					for (int i = 0; i < spectrum.length; i ++) {
						for (int j = 0; j < spectrum[i].length; j ++) { // compute the error between it and the actual spectrum
							if (fitSpectrum[i][j] + noiseVar[j] > 0) // skipping places where we expect 0 and got 0
								err += Math.log(fitSpectrum[i][j] + noiseVar[j])/2 +
										Math.pow(fitSpectrum[i][j] - spectrum[i][j], 2)/
												(2*(fitSpectrum[i][j] + noiseVar[j]));
							else if (spectrum[i][j] != 0) // but throwing infinity if we expect 0 and got not 0
								err = Double.POSITIVE_INFINITY;
						}
					}
					assert Double.isFinite(err) : err;
					
					double penalty = 0; // negative log of prior (ignoring global normalization)
					for (int i = 0; i < spectrum.length; i ++) {
						for (int j = 0; j < spectrum[i].length; j ++) {
							if (teoSpectrum[i][j] > 1e-300)
								penalty -= 1.0*teoSpectrum[i][j]*
										(1 - Math.log(teoSpectrum[i][j]/spectrumScale))/spectrumScale; // encourage entropy
							else
								penalty -= 0;
						}
					}
					final double dt = timeAxis[1] - timeAxis[0];
					for (int k = 0; k < 4; k ++) {
						double curveScale = 8*PARAM_SCALES[k]/Math.pow(MAX_T - MIN_T, 2); // discourage egregious curvature
						if (k == 0)
							curveScale = 8*spectrumScale/Math.pow(MAX_T - MIN_T, 2);
						for (int i = 1; i < timeAxis.length - 1; i ++) {
							penalty += 0.01*Math.pow((params[k][i-1] - 2*params[k][i] + params[k][i+1])/(dt*dt)/curveScale, 2);
						}
					}
					
//					System.out.println(err+penalty);
					return err + penalty;
				})
		).getPoint();
//				}, initialGuess, lowerBounds, upperBounds, 1e-14);
		
		long endTime = System.currentTimeMillis();
		if (logger != null)
			logger.info(String.format(Locale.US, "completed in %.2f minutes.",
					(endTime - startTime)/60000.));
		
		double[][] params = new double[4][timeAxis.length];
		for (int k = 0; k < 4; k ++) // first unpack the state vector
			for (int i = 0; i < timeAxis.length; i ++)
				params[k][i] = opt[i+k*timeAxis.length];
		this.neutronYield   = params[0];
		this.ionTemperature = params[1];
		this.flowVelocity   = params[2];
		this.arealDensity   = params[3];
		
		this.fitNeutronSpectrum = generateSpectrum( // and then interpret it
				neutronYield, ionTemperature, flowVelocity, arealDensity, energyBins, timeBins);
		this.fitDeuteronSpectrum = this.response(energyBins, timeBins, fitNeutronSpectrum, false);
		
		final double dt = (timeBins[1] - timeBins[0]);
		for (int i = 0; i < timeAxis.length; i ++) {
			if (neutronYield[i]*1e15*efficiency(14.1e6)*dt < MIN_STATISTICS) { // now remove any intensive values from times with insufficient statistics
				ionTemperature[i] = Double.NaN;
				flowVelocity[i] = Double.NaN;
				arealDensity[i] = Double.NaN;
			}
		}
		
		double[] dTidt = NumericalMethods.derivative(timeAxis, ionTemperature); // now we can freely analyze the resulting profiles
		double[] dρRdt = NumericalMethods.derivative(timeAxis, arealDensity);
		double[] dvidt = NumericalMethods.derivative(timeAxis, flowVelocity);
		int bangTime = NumericalMethods.argmax(neutronYield);
		int compressTime = NumericalMethods.argmax(arealDensity);
		int maxPRRamp = NumericalMethods.argmax(dρRdt);
		double[] moments = new double[5];
		for (int k = 0; k < moments.length; k ++)
			moments[k] = NumericalMethods.moment(k, timeBins, neutronYield);
		
		if (bangTime == -1 || compressTime == -1 || maxPRRamp == -1) {
			if (logger != null)
				logger.warning("Insufficient statistics to analyze.");
			return null;
		}
		else {
			if (logger != null) {
				logger.info(String.format("Bang time:         %8.3f ns", timeAxis[bangTime]));
				logger.info(String.format("Peak compression:  %8.3f ns", timeAxis[compressTime]));
				logger.info(String.format("            = BT + %8.3f ps", (timeAxis[compressTime] - timeAxis[bangTime])/1e-3));
				logger.info(String.format("Max ρR ramp:       %8.3f ps", timeAxis[maxPRRamp]));
				logger.info(String.format("            = BT + %8.3f ps", (timeAxis[maxPRRamp] - timeAxis[bangTime])/1e-3));
				logger.info(String.format("Ti at BT:          %8.3f keV", ionTemperature[bangTime]));
				logger.info(String.format("vi at BT:          %8.3f μm/ns", flowVelocity[bangTime]));
				logger.info(String.format("dTi/dt at BT:      %8.3f keV/ns", dTidt[bangTime]));
				logger.info(String.format("dρR/dt at BT:      %8.3f g/cm^2/ns", dρRdt[bangTime]));
				logger.info(String.format("dvi/dt at BT:      %8.3f μm/ns^2", dvidt[bangTime]));
				logger.info(String.format("Peak ρR:           %8.3f g/cm^2", arealDensity[compressTime]));
				logger.info(String.format("Total yield (μ0):  %8.3g", moments[0]));
				logger.info(String.format("Burn mean (μ1):    %8.3g ns", moments[1]));
				logger.info(String.format("Burn width (μ2):   %8.3g ns", Math.sqrt(moments[2])*2.355));
				logger.info(String.format("Burn skewness (μ3):%8.3g", moments[3]));
				logger.info(String.format("Burn kurtosis (μ4):%8.3g", moments[4]));
			}
			return new double[] {
					(endTime - startTime)/1000., timeAxis[bangTime], timeAxis[compressTime], timeAxis[maxPRRamp],
					ionTemperature[bangTime], arealDensity[bangTime], flowVelocity[bangTime], dTidt[bangTime],
					dρRdt[bangTime], dvidt[bangTime], arealDensity[compressTime], moments[0],
					moments[1], Math.sqrt(moments[2])*2.355, moments[3], moments[4]
			};
		}
	}
	
	/**
	 * simulate a single random neutron emitted from TCC at the given energy and determine the
	 * position and time at which its child ion crosses the focal plane.
	 * @param energy initial energy of released neutron [eV].
	 * @param time initial time of released neutron [s].
	 * @return { signed hypot(x,z), t } [m, s].
	 */
	public double[] simulate(double energy, double time) {
		double[] rCollision = chooseCollisionPosition();
		
		double[] rAperture = chooseAperturePosition();
		
		double[] vFinal = computeFinalVelocity(energy, rCollision, rAperture);
		
		double[] rFocal = computeFocusedPosition(rCollision, vFinal, time);
		
		return new double[] { rFocal[x]/Math.cos(focalPlaneAngle), rFocal[3] };
	}
	
	/**
	 * estimate the original time and energy of this ion's neutron without looking at its actual
	 * time and energy, by guessing its energy and accounting for travel time.
	 * @param position the position where it hits the focal plane [m]
	 * @param time the time at which it hits the focal plane [s]
	 * @return { energy, time } [J, s].
	 */
	public double[] backCalculate(double position, double time) {
		double focusingDistance = position*Math.sin(focalPlaneAngle);
		double E = energyVsPosition.evaluate(position); // [J]
		double v = Math.sqrt(2*E/ion.mass);
		double d0 = (E - cosyK0)/cosyK0;
		double lf = cosyPolynomial(4, new double[] {0, 0, 0, 0, 0, d0}); // re-use the COSY mapping to estimate the time of flight from energy
		double timeCorrection = cosyT0 + lf*cosyT1 + focusingDistance/v;
		return new double[] { E/energyFactor, time - timeCorrection };
	}
	
	/**
	 * choose a random location in the foil for the neutron to collide.
	 * @return { x, y, z } [m]
	 */
	private double[] chooseCollisionPosition() {
		double θF = Math.acos(1 - 2*Math.random()*probHitsFoil);
		double rF = foilDistance*Math.tan(θF); // NOTE: original code assumes uniform distribution within foil; I account for nonzero solid angle subtended at TCC.
		double φF = Math.random()*2*Math.PI;
		double zF = foilDistance + (2*Math.random()-1)*foilThickness/2; // assume foil is thin, so every z coordinate is equally likely
		return new double[] { rF*Math.cos(φF), rF*Math.sin(φF), zF };
	}
	
	/**
	 * choose a random location in the aperture plane for the deuteron to pass through.
	 * @return { x, y, z } [m]
	 */
	private double[] chooseAperturePosition() {
		double xA = (2*Math.random()-1)*apertureWidth/2; // assume aperture is far away, so every point in it is equally likely to be hit
		double yA = (2*Math.random()-1)*apertureHeight/2;
		double zA = apertureDistance;
		return new double[] { xA, yA, zA };
	}
	
	/**
	 * compute the velocity with which the deuteron passes through the aperture.
	 * @param vInitial {vx,vy,vz} of the neutron as it enters the foil [m/s].
	 * @param A ratio of charged particle mass to neutron mass.
	 * @param rFoil {x,y,z} of the point at which the neutron strikes the deuteron [m].
	 * @param rAperture {x,y,z} of the point where the deuteron passes through the aperture [m].
	 * @return { vx, vy, vz } [m/s]
	 */
	private double[] computeFinalVelocity(
			double energy, double[] rFoil, double[] rAperture) {
		double E0 = -Particle.E.charge*energy; // convert energy from [eV] to [J]
		
		double[] nHat = { // get the unit vector in the direction of the neutron
				rFoil[x], rFoil[y], rFoil[z] };
		double norm = Math.sqrt(sqr(nHat));
		for (int i = 0; i < 3; i ++)
			nHat[i] /= norm;
		
		double[] dHat = { // and the unit vector in the direction of the ion
				rAperture[0] - rFoil[0], rAperture[1] - rFoil[1], rAperture[2] - rFoil[2] };
		norm = Math.sqrt(sqr(dHat));
		for (int i = 0; i < 3; i ++)
			dHat[i] /= norm;
		
		double cosθ = nHat[x]*dHat[x] + nHat[y]*dHat[y] + nHat[z]*dHat[z];
		double E1 = energyFactor*E0*cosθ*cosθ; // assume elastic collision between neutron and ion
		double distance = (foilDistance + foilThickness/2 - rFoil[2])/dHat[2];
		E1 = energyVsDistance.evaluate(distanceVsEnergy.evaluate(E1) - distance); // lose some energy by dragging through the foil
		
//		System.out.print(E1/Particle.P.charge/1e6+", ");
		
		if (E1 < cosyKmin || E1 > cosyKmax) {
			return new double[] { Double.NaN, Double.NaN, Double.NaN }; // some won't make it through the "energy aperture"
		}
		else {
			double v = Math.sqrt(2*E1/ion.mass); // get the resulting velocity otherwise
			return new double[] { v*dHat[x], v*dHat[y], v*dHat[z] };
		}
	}
	
	/**
	 * evaluate the mapping provided by COSY to obtain the time and position and velocity at
	 * which the ion passes the back reference plane.
	 * @param rFoil {x,y,z} of the point at which the neutron strikes the deuteron [m].
	 * @param vInit {vx,vy,vz} of the deuteron as it exits the foil [m/s]
	 * @param tNeutron the time at which the neutron was released [s].
	 * @return { x, y, z, t } at which it strikes the focal plane [m, m, s]
	 */
	private double[] computeFocusedPosition(double[] rFoil, double[] vInit, double tNeutron) {
		if (Double.isNaN(vInit[0]))
			return new double[] { Double.NaN, Double.NaN, Double.NaN, Double.NaN };
		
		double x0 = rFoil[x], y0 = rFoil[y]; // COSY takes spatial coordinates in [m] (assume foil is thin so we can ignore rFoil[2])
		double a0 = vInit[x]/cosyV0, b0 = vInit[y]/cosyV0; // angular coordinates in [rad] (more or less)
		double t0 = 0; // assume time it takes neutron to hit foil is negligible
		double K0 = 1/2.*ion.mass*sqr(vInit); // for the 'd' coordinate, we must convert velocity to energy
		double d0 = (K0 - cosyK0)/cosyK0; // and then compare that to an expected reference energy
		double[] input = { x0, a0, y0, b0, t0, d0 };
		double[] output = new double[5];
		for (int i = 0; i < output.length; i ++) // the polynomial is pretty simple to compute
			output[i] = cosyPolynomial(i, input);
		double[] rPlane = { output[0], output[2], 0 }; // [m]
		double[] vFinal = { output[1]*cosyV0, output[3]*cosyV0, vInit[z] }; // [m/s]
		double tPlane = tNeutron + cosyT0 + output[4]*cosyT1;
		
		double tFocusing = rPlane[x] / (vFinal[z]/Math.tan(focalPlaneAngle) - vFinal[x]); // finally, account for the additional t that it takes to strike the focal plane
		double[] rFocused = {
				rPlane[x] + tFocusing*vFinal[x],
				rPlane[y] + tFocusing*vFinal[y],
				rPlane[z] + tFocusing*vFinal[z] };
		
		return new double[] { rFocused[x], rFocused[y], rFocused[z], tPlane + tFocusing };
	}
	
	/**
	 * evaluate a single polynomial from that table of coefficients.
	 * @param i the index of the parameter to compute
	 * @param input the initial values to use
	 * @return the final value of parameter i
	 */
	private double cosyPolynomial(int i, double[] input) {
		double output = 0;
		for (int j = 0; j < cosyCoefficients.length; j ++) {
			double term = cosyCoefficients[j][i];
			for (int k = 0; k < input.length; k ++) {
				term *= Math.pow(input[k], cosyExponents[j][k]);
			}
			output += term;
		}
		return output;
	}
	/**
	 * generate a time-resolved spectrum based on some parameters that vary with time.
	 * @param Yn the neutron yield rate [10^15/ns]
	 * @param Ti the ion temperature [keV]
	 * @param vi the bulk flow rate parallel to the line of sight [μm/ns]
	 * @param ρR the areal density of fuel and shell surrounding the hot spot [g/cm^2]
	 * @param t the edges of the time bins [ns]
	 * @param eBins the edges of the energy bins [MeV]
	 * @return the theoretical number of particles in each energy bin, ignoring stochastity.
	 */
	public static double[][] generateSpectrum(double[] Yn, double Ti[], double vi[], double ρR[],
			double[] eBins, double[] tBins) {
		double[][] spectrum = new double[eBins.length-1][tBins.length-1];
		for (int j = 0; j < spectrum[0].length; j ++) {
			double dt = (tBins[j+1] - tBins[j]); // bin width [ns]
			double[] timeSlice = generateSpectrum(Yn[j]*dt, Ti[j], vi[j], ρR[j], eBins);
			for (int i = 0; i < spectrum.length; i ++)
				spectrum[i][j] = timeSlice[i];
		}
		return spectrum;
	}
	
	/**
	 * generate a time-averaged spectrum based on some parameters that are taken to be constant.
	 * @param Yn the total neutron yield [10^15]
	 * @param Ti the ion temperature [keV]
	 * @param vi the bulk flow rate parallel to the line of sight [μm/ns]
	 * @param ρR the areal density of fuel and shell surrounding the hot spot [g/cm^2]
	 * @param eBins the edges of the energy bins [MeV]
	 * @return the theoretical number of particles in each energy bin, ignoring stochastity.
	 */
	public static double[] generateSpectrum(double Yn, double Ti, double vi, double ρR,
			double[] eBins) {
		double Ps = 1 - Math.exp(-0.212*ρR); // probability of any scattering or absorption (1.80b/5u) []
		double p = .0590*ρR; // DS spectrum parameter (.493b/MeV)/5u [MeV^-1]
		double λ = 0.313; // DS spectrum parameter (1/3.19MeV) [MeV^-1]
		double μ = 14.1 + .54e-3*vi; // primary peak (see paper), [MeV]
		double s2 = 0.403e-3*μ*Ti; // primary variance (see paper) [MeV^2]
		double[] pdf = new double[eBins.length]; // probability distribution at edges
		double[] counts = new double[eBins.length-1]; // spectrum counts in bins
		for (int i = eBins.length-1; i >= 0; i --) {
			if (Ti > 0) {
				double E = eBins[i]; // x coordinate [MeV]
				double primaryDistribution = (1 - Ps)/Math.sqrt(2*Math.PI*s2)*
						Math.exp(-Math.pow(E - μ, 2)/(2*s2)); // primary distribution [MeV^-1]
				double dsDistribution = p/2*Math.exp(λ*(E - μ + s2*λ/2))*
						(1 + NumericalMethods.erf((μ - E + s2*λ)/Math.sqrt(2*s2))); // secondary distribution [MeV^-1]
				pdf[i] = primaryDistribution + dsDistribution; // combine the two distributions
			}
			else
				pdf[i] = 0;
			if (i < eBins.length-1) // if we have the pdf at both edges of this bin
				counts[i] = Yn*1e15*(pdf[i] + pdf[i+1])/2*(eBins[i+1] - eBins[i]); // fill it trapezoidally []
		}
		return counts;
	}
	
	public double[] getTimeBins() {
		return this.timeBins;
	}
	
	public double[] getEnergyBins() {
		return this.energyBins;
	}
	
	public double[][] getCorrectedSpectrum() {
		return this.deuteronSpectrum;
	}
	
	public double[][] getInferredSpectrum() {
		return this.fitNeutronSpectrum;
	}
	
	public double[][] getFittedSpectrum() {
		return this.fitDeuteronSpectrum;
	}
	
	public double[] getTimeAxis() {
		return this.timeAxis;
	}
	
	public double[] getIonTemperature() {
		return this.ionTemperature;
	}
	
	public double[] getFlowVelocity() {
		return this.flowVelocity;
	}
	
	public double[] getArealDensity() {
		return this.arealDensity;
	}
	
	public double[] getNeutronYield() {
		return this.neutronYield;
	}
	
	/**
	 * convert a time-cumulative spectrum, as read from an input file, into an array of counts
	 * densities. also, edit energies to make it bin edges. this will remove the last row of
	 * energies. given that the last row is over 29 MeV, I think that's fine.
	 * @param timeAxis the time bin boundaries [ns]
	 * @param energyAxis the energy bin midpoints [MeV]
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
	 * square a vector, because this takes so long to write out every time.
	 * @param v
	 * @return
	 */
	private static double sqr(double[] v) {
		double s = 0;
		for (double x: v)
			s += Math.pow(x, 2);
		return s;
	}
	
	
	/**
	 * @param args
	 * @throws IOException 
	 * @throws NumberFormatException 
	 */
	public static void main(String[] args) throws NumberFormatException, IOException {
//		MRSt sim = new MRSt(
//				Particle.D, 3.0e-3, 0.3e-3, 80e-6,
//				CSV.read(new File("data/stopping_power_deuterons.csv"), ','),
//				6e0, 4.0e-3, 20.0e-3, 10.7e6, 14.2e6, 12.45e6,
//				CSV.readCosyCoefficients(new File("data/MRSt_IRF_FP tilted.txt"), 3),
////				CSV.readCosyCoefficients(new File("data/MRSt_IRF_FP not tilted.txt"), 3),
//				CSV.readCosyExponents(new File("data/MRSt_IRF_FP tilted.txt"), 3), 70.3,
//				100, null);
//		
//		double[][] spectrum = CSV.read(new File("data/nsp_150327_16p26.txt"), '\t');
//		double[] timeAxis = CSV.readColumn(new File("data/nsp_150327_16p26_time.txt"));
//		double[] energyAxis = CSV.readColumn(new File("data/Energy bins.txt"));
//		spectrum = interpretSpectrumFile(timeAxis, energyAxis, spectrum);
		
//		for (int i = 0; i < 216; i ++) {
//			double e0 = 12+Math.random()*5;
//			System.out.print("["+e0+", ");
//			double[] xt = sim.simulate(e0*1e6, 0);
//			System.out.print(xt[0]/1e-2+", ");
//			double[] et = sim.backCalculate(xt[0], xt[1]);
//			double e = et[0]/(-Particle.E.charge)/1e6, t0 = et[1]/1e-9;
//			System.out.println(e+", "+t0+"],");
//		}
		
		int n = 40;
		double[] E = new double[n+1];
		for (int i = 0; i < n+1; i ++)
			E[i] = 12 + 4./n*i;
		System.out.println(Arrays.toString(E));
		System.out.println(Arrays.toString(generateSpectrum(10, 2, 100, 0.5, E)));
		
//		sim.respond(energyAxis, timeAxis, spectrum);
	}
}
