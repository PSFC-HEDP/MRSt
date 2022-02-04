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
package physics;

import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.MaxIter;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.MultivariateOptimizer;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.MultiDirectionalSimplex;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.PowellOptimizer;
import util.COSYMapping;
import util.NumericalMethods;
import util.NumericalMethods.Quantity;
import util.Optimization;

import java.io.IOException;
import java.util.Arrays;
import java.util.Locale;
import java.util.Random;
import java.util.function.Function;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * the class where all the math is.
 * 
 * @author Justin Kunimune
 */
public class Analysis {
	
	public static final String[] HEADERS = {
		  "Yield factor", "Temperature factor", "Down-scatter factor", "Velocity shift (μm/ns)",
		  "Computation time (s)", "Total yield (10^15)", "Bang time (ns)",
		  "Burn width (ns)", "Burn skewness", "Burn kurtosis", "Stagnation - BT (ns)",
		  "Burn-average Ti (keV)", "Peak Ti (keV)",
		  "Ti at BT (keV)", "dTi/dt at BT (keV/ns)", "d^2Ti/dt^2 at BT (keV/ns^2)",
		  "Burn-average \u03C1R (g/cm^2)", "\u03C1R at stagnation (g/cm^2)",
		  "\u03C1R at BT (g/cm^2)", "d\u03C1R/dt at BT (g/cm^2/ns)", "d^2V/dt^2/V at BT (1/ns^2)",
		}; // the names, units, and order of time-dependent burn parameters
	public static final String[] HEADERS_WITH_ERRORS = appendErrorsToHeader();

	public static final Random MC_RANDOM = new Random(1);
	public static final Random NOISE_RANDOM = new Random(0);

	private static final double MIN_E = 12, MAX_E = 16; // histogram bounds [MeV]
	private static final int BUFFER = 4; // empty pixels to include simulate on each side [ns]
	private static final double E_BIN = .05, T_BIN = 30e-3; // bin sizes [MeV], [ns]

//	private static final double SUBSTRATE_THICKNESS = 100; // [μm]
//	private static final double PHOTOCATHODE_THICKNESS = .1; // [μm]
//	private static final double PDDT_BIAS = 1e3; // [V]
//	private static final double MESH_LENGTH = 1e-3; // [m]
//	private static final double DRIFT_LENGTH = 1e0; // [m]
//	private static final double TIME_DILATION = 20;
//	private static final double MCT_POROSITY = .70;
//	private static final double MCT_GAIN = 1e4;
//
//	private static final double TRANSFER_FUNC_ERROR = 0.00; // the error in the transfer function

	private final IonOptics ionOptics; // the ion optic system
	private final Detector detector; // the detector system

	private final double precision; // factor by which to ease the convergence conditions

	private final double[] energyBins; // endpoints of E bins for inferred spectrum [MeV]
	private double[] timeBins; // endpoints of time bins for inferred spectrum [ns]
	private double[][] deuteronSpectrum; // time-corrected deuteron counts
	private double[][] fitNeutronSpectrum; // backward-fit neutron counts
	private double[][] fitDeuteronSpectrum; // backward-fit deuteron counts (this should be similar to deuteronSpectrum)
	
	private final double[] energyAxis; // centers of energy bins [ns]
	private double timeStep;
	private double[] timeAxis; // centers of time bins [ns]
	private Quantity[] neutronYield; // 1e15/ns
	private Quantity[] ionTemperature; // keV
	private Quantity[] arealDensity; // g/cm^2
	private double[][] covarianceMatrix; // and covariances that go with all of these
	
	private final Logger logger; // for logging
	
	/**
	 * perform some preliminary calculations for the provided configuration.
	 * @param foilDistance the distance from TCC to the foil [m]
	 * @param foilWidth the total width of the foil [m]
	 * @param foilHeight the total hite of the foil [m]
	 * @param foilThickness the thickness of the foil [m]
	 * @param apertureDistance the distance from TCC to the aperture [m]
	 * @param apertureWidth the width of the aperture [m]
	 * @param apertureHeight the hite of the aperture [m]
	 * @param cosyMapping the COSY matrix
	 * @param focalTilt angle of the focal plane (0 means untilted) [deg]
	 * @param precision a number that modifies how hard it tries to fully converge
	 * @throws IOException if one or more of the stopping power tables cannot be
	 *                     accessd for any reason
	 */
	public Analysis(
			double foilDistance, double foilWidth, double foilHeight, double foilThickness,
			double apertureDistance, double apertureWidth, double apertureHeight,
			COSYMapping cosyMapping, double focalTilt, boolean reuseMatrix,
			double precision, Logger logger) throws IOException {

		this.ionOptics = new IonOptics(
				foilDistance, foilWidth, foilHeight, foilThickness,
				apertureDistance, apertureWidth, apertureHeight,
				MIN_E, MAX_E, cosyMapping, focalTilt, reuseMatrix);
//		this.detector = new PulseDilationDriftTube(
//				ion, SUBSTRATE_THICKNESS, PHOTOCATHODE_THICKNESS,
//				focalTilt, PDDT_BIAS, MESH_LENGTH, DRIFT_LENGTH,
//				TIME_DILATION, MCT_POROSITY, MCT_GAIN, 100);
		this.detector = new StreakCameraArray(
			  2.5e-2, 400e-6,
			  (focalTilt == 0) ? 11.5e-9 : 4.5e-9,
			  2.4/Math.cos(Math.toRadians(focalTilt)) * 51,
			  81/(25e-6*25e-6),
			  81/(25e-6*25e-6),
			  (focalTilt == 0) ? new double[] {-0.75e-2} : new double[] {-5e-2, 0, 5e-2},
			  ionOptics);
//		this.detector = new PerfectDetector();

		this.precision = precision;

		this.energyBins = new double[(int) ((MAX_E - MIN_E)/E_BIN + 1)];
		for (int i = 0; i < energyBins.length; i ++)
			this.energyBins[i] = (MIN_E + i*(MAX_E - MIN_E)/(energyBins.length-1));
		
		this.energyAxis = new double[energyBins.length-1];
		for (int i = 0; i < energyBins.length-1; i ++)
			this.energyAxis[i] = (this.energyBins[i] + this.energyBins[i+1])/2;

		this.logger = logger;
	}
	
	/**
	 * establish the time bins and related quantities
	 * @param spectrumTimeBins the time bins of the input spectrum, to use as
	 *                         a reference
	 */
	private void instantiateTimeAxis(double[] spectrumTimeBins) {
		double minT = spectrumTimeBins[0] - BUFFER*T_BIN;
		double maxT = spectrumTimeBins[spectrumTimeBins.length-1] + BUFFER*T_BIN;
		this.timeBins = new double[(int) ((maxT - minT)/T_BIN + 1)];
		for (int i = 0; i < timeBins.length; i ++)
			this.timeBins[i] = minT + i*(maxT - minT)/(timeBins.length-1);
		
//		double[][] actual;
//		try {
//			actual = CSV.read(new File("input/trajectories og with falling temp.csv"), ',', 1);
//		} catch (NumberFormatException e) {
//			actual = null;
//		} catch (IOException e) {
//			actual = null;
//		}
//		this.timeBins = new double[actual.length+1];
//		for (int i = 0; i < timeBins.length-1; i ++)
//			this.timeBins[i] = actual[i][0] - (actual[1][0] - actual[0][0])/2.;
//		this.timeBins[timeBins.length-1] = actual[actual.length-1][0] + (actual[1][0] - actual[0][0])/2.;
		
		this.timeStep = timeBins[1] - timeBins[0];
		this.timeAxis = new double[timeBins.length-1];
		for (int i = 0; i < timeBins.length-1; i ++)
			this.timeAxis[i] = (this.timeBins[i] + this.timeBins[i+1])/2;
	}

	/**
	 * compute the dispersion at the detector plane
	 * @return dE/dl [keV/mm]
	 */
	public double computeDispersion() {
		return this.ionOptics.computeDispersion();
	}

	/**
	 * compute the time skew at the detector plane
	 * @return dt/dE [ps/keV]
	 */
	public double computeTimeSkew() {
		return this.ionOptics.computeTimeSkew();
	}

	/**
	 * compute the time and energy resolution (FWHM) for particles at a given energy.
	 * @param referenceEnergy the energy of the central ray [MeV]
	 * @return {energy resolution [keV], time resolution [ps]}
	 */
	public double[] computeResolution(double referenceEnergy) {
		return this.ionOptics.computeResolution(referenceEnergy);
	}

	public double efficiency() {
		return this.ionOptics.efficiency(14);
	}

	public double gain() {
		return this.ionOptics.efficiency(14)
				*this.detector.efficiency(14)
				*this.detector.gain();
	}

	/**
	 * what does the detector see if the given neutron spectrum goes thru the
	 * ion optics
	 */
	public double[][] response(double[] energyBins, double[] timeBins,
							   double[][] spectrum,
							   boolean stochastic, boolean actual) {
//		System.out.println("responding...");
//		System.out.println(NumericalMethods.sum(spectrum));
//		System.out.println(NumericalMethods.sum(this.ionOptics.response(energyBins, timeBins, spectrum, stochastic, actual)));
//		System.out.println(NumericalMethods.sum(this.detector.response(energyBins, timeBins, ionOptics.response(energyBins, timeBins, spectrum, false, actual), false)));
		return this.detector.response(
			  energyBins, timeBins,
			  this.ionOptics.response(
			  	  energyBins, timeBins, spectrum, stochastic, actual),
			  stochastic);
	}
	
	/**
	 * compute the total response to an implosion with the given neutron spectrum
	 * and save and analyze it.
	 * @param energies the energies that describe the rows of counts [MeV]
	 * @param times the times that describe the columns of counts [ns]
	 * @param spectrum the time- and energy- resolved neutron spectrum in number
	 *                 of neutrons. each row corresponds to one element of
	 *                 energies, and each column one element of times. [#/MeV/ns]
	 * @param errorBars whether to bother computing error bars
	 * @return {computation time, 0, Yn, err, BT, err, BW, err, skewness, err,
	 *          kurtosis, err, peak compression, err, rho R (BT), err,
	 *          \< rho R \>, err, d(rho R)/dt, err, \<Ti\>, err, dTi/dt, err,
	 *          \<vi\>, err, dvi/dt, err}
	 */
	public double[] respondAndAnalyze(double[] energies, double[] times, double[][] spectrum, ErrorMode errorBars) {
		if (spectrum.length != energies.length-1 || spectrum[0].length != times.length-1)
			throw new IllegalArgumentException("These dimensions don't make any sense.");

		this.instantiateTimeAxis(times); // first, create the detector time bins

		spectrum = NumericalMethods.downsample(
				times, energies, spectrum, timeBins, energyBins); // first rebin the spectrum

		logger.info("beginning Monte Carlo computation.");
		long startTime = System.currentTimeMillis();
		this.deuteronSpectrum = this.response(
				energyBins, timeBins, spectrum, true, true);
		long endTime = System.currentTimeMillis();
		logger.info(String.format(Locale.US, "completed in %.2f minutes.",
			                      (endTime - startTime)/60000.));

		Quantity[] values = analyze(deuteronSpectrum, errorBars);

		if (values != null) {
			double[] output = new double[2*values.length];
			for (int i = 0; i < values.length; i++) {
				output[2*i + 0] = values[i].value;
				output[2*i + 1] = Math.sqrt(values[i].variance(covarianceMatrix));
			}
			return output;
		}
		else {
			return null;
		}
	}

	/**
	 * take the time-corrected spectrum and use it to compute and store time-resolved values
	 * for ion temperature, areal density, and yield.
	 * @param spectrum the time-corrected spectrum we want to understand
	 * @param errorMode should error bars be computed (it's rather intensive)?
	 * @return {computation time, Yn, BT, BW, skewness, kurtosis, peak compression, rho R (BT),
	 * \< rho R \>, d(rho R)/dt, \<Ti\>, dTi/dt, \<vi\>, dvi/dt} or null if it can't even
	 */
	private Quantity[] analyze(double[][] spectrum, ErrorMode errorMode) {
		if (spectrum.length != energyBins.length-1 || spectrum[0].length != timeBins.length-1)
			throw new IllegalArgumentException("These dimensions are wrong.");

		if (NumericalMethods.max(spectrum) == 0) {
			logger.log(Level.SEVERE, "There were no deuterons detected.");
			return null;
		}
		else {
			logger.log(Level.INFO, String.format("There were %.1g deuterons detected.", NumericalMethods.sum(spectrum)));
		}

		logger.info("beginning fit process.");
		long startTime = System.currentTimeMillis();

		for (double[] doubles: spectrum)
			for (double value: doubles)
				if (Double.isNaN(value))
					throw new IllegalArgumentException("well fuck you too");

		final int N = 3;
		final int M = timeAxis.length;

		// start with the simplest possible reconstruction
		double[] naiveNeutronYield = new double[spectrum[0].length];
		for (int i = 0; i < spectrum.length; i ++) {
			double background = this.detector.background(energyAxis[i], energyBins, timeBins);
			for (int j = 0; j < spectrum[i].length; j ++) {
				naiveNeutronYield[j] += Math.max(
					  0, (spectrum[i][j] - background)/this.gain())/(timeStep*1e15); // the spectrum scaled by the efficiency of the system
			}
		}
		//		double[][] D = new double[energyBins.length-1][timeBins.length-1];
		//		for (int i = 0; i < energyBins.length-1; i ++)
		//			for (int j = 0; j < timeBins.length-1; j ++)
		//				D[i][j] = Math.max(1, spectrum[i][j]);
		//		double[][] gelf = Optimization.optimizeGelfgat(spectrum, D, this.rongTransferMatrix, 1e5);

		double[] yieldScale = new double[M];
		double[] temperatureScale = new double[M];
		double[] densityScale = new double[M];
		double[] lowerBound = new double[M];
		double[] upperBound = new double[M];
		for (int j = 0; j < M; j++) {
			yieldScale[j] = NumericalMethods.max(naiveNeutronYield)/6.;
			temperatureScale[j] = 10;
			densityScale[j] = 1;
			lowerBound[j] = 0;
			upperBound[j] = Double.POSITIVE_INFINITY;
		}

		// put together the initial gess for the proper fitting
		double[] ionTemperature = new double[timeAxis.length]; // keV
		double[] electronTemperature = new double[timeAxis.length]; // keV
		double[] bulkFlowVelocity = new double[timeAxis.length]; // μm/ns
		double[] arealDensity = new double[timeAxis.length]; // g/cm^2
		for (int j = 0; j < timeAxis.length; j ++) {
			ionTemperature[j] = 4;
			electronTemperature[j] = 4;
			bulkFlowVelocity[j] = 0;
			arealDensity[j] = 200e-3;
		}

		// then do the second simplest reconstruction for the burn history
		double[] finalIonTemperature = ionTemperature;
		double[] finalArealDensity = arealDensity;
		double[] neutronYield = optimize(
			  (double[] x) -> logPosterior(
					spectrum,
					x,
					finalIonTemperature,
					electronTemperature,
					bulkFlowVelocity,
					finalArealDensity,
					null, true, 1),
			  naiveNeutronYield, yieldScale, lowerBound, upperBound); // 10^15/ns

		// determine the bounds of the region that's worth optimizing
		double[] statistics = new double[M];
		for (int i = 0; i < M; i ++)
			statistics[i] = neutronYield[i]*1e15*this.efficiency()*timeStep;
		final int peak = NumericalMethods.argmax(statistics);
		final int left = Math.min(peak - 2, NumericalMethods.firstAbove(1, statistics));
		final int rite = Math.max(peak + 3, NumericalMethods.lastAbove(1, statistics) + 1);
		boolean[] active = new boolean[timeAxis.length];
		for (int j = 0; j < timeAxis.length; j ++)
			active[j] = j >= left && j < rite;

		//		for (double smoothing = 1e4; smoothing > .9; smoothing /= 10) {
		final double finalSmoothing = 1;//smoothing;

		double lastValue;
		double thisValue = logPosterior(
			  spectrum, neutronYield,
			  ionTemperature, electronTemperature,
			  bulkFlowVelocity, arealDensity,
			  active, false, finalSmoothing);

		do {
			final double[] neutronYieldInitialGess = neutronYield;
			final double[] ionTemperatureInitialGess = ionTemperature;
			final double[] arealDensityInitialGess = arealDensity;

			logger.log(Level.FINE, "Fitting ion temperature...");

			final double[] ionTemperatureNewGess = optimize(
				  (double[] x) -> logPosterior(
						spectrum,
						neutronYieldInitialGess,
						x,
						electronTemperature,
						bulkFlowVelocity,
						arealDensityInitialGess,
						active, false, finalSmoothing),
				  ionTemperatureInitialGess,
				  temperatureScale, lowerBound, upperBound, active);

			logger.log(Level.FINE, "Fitting ρR trajectory...");

			final double[] arealDensityNewGess = optimize(
				  (double[] x) -> logPosterior(
					  spectrum,
					  neutronYieldInitialGess,
					  ionTemperatureNewGess,
					  electronTemperature,
					  bulkFlowVelocity,
					  x,
					  active, false, finalSmoothing),
				  arealDensityInitialGess,
				  densityScale, lowerBound, upperBound, active);

			logger.log(Level.FINE, "Fitting burn history...");

			final double[] neutronYieldNewGess = optimize(
				  (double[] x) -> logPosterior(
						spectrum,
						x,
						ionTemperatureNewGess,
						electronTemperature,
						bulkFlowVelocity,
						arealDensityNewGess,
						active, false, finalSmoothing),
				  neutronYieldInitialGess,
				  yieldScale, lowerBound, upperBound, active);

			ionTemperature = ionTemperatureNewGess;
			arealDensity = arealDensityNewGess;
			neutronYield = neutronYieldNewGess;

			lastValue = thisValue;
			thisValue = logPosterior(
				  spectrum, neutronYield,
				  ionTemperature, electronTemperature,
				  bulkFlowVelocity, arealDensity,
				  active, false, finalSmoothing);
		} while (lastValue - thisValue > this.precision);
//		}

//		try {
//			double ultimatePosterior = this.logPosterior(spectrum, neutronYield, ionTemperature, electronTemperature, bulkFlowVelocity, arealDensity, left, rite, false, finalSmoothing);
//			double ultimateBayesor = this.logPosterior(spectrum, neutronYield, ionTemperature, electronTemperature, bulkFlowVelocity, arealDensity, left, rite, false, 0);
//			double[][] trueTrajectories;
//			trueTrajectories = CSV.read(new File("input/trajectories og with falling temp.csv"), ',', 1);
//			double[] time = new double[trueTrajectories.length];
//			double[] ρR = new double[trueTrajectories.length];
//			double[] Yn = new double[trueTrajectories.length];
//			double[] Ti = new double[trueTrajectories.length];
//			for (int i = 0; i < trueTrajectories.length; i++) {
//				time[i] = trueTrajectories[i][0];
//				Yn[i] = trueTrajectories[i][1]*.1*1e6/1e-6/(14e6*1.6e-19)/1e15*1e-9; // convert from 0.1MJ/μs to 1e15n/ns
//				Ti[i] = trueTrajectories[i][4];
//				ρR[i] = trueTrajectories[i][3];
//			}
//			for (int j = NumericalMethods.argmax(neutronYield); j <= NumericalMethods.argmax(neutronYield)+1; j ++) {
//				if (timeAxis[j] < time[0]) {
////					neutronYield[j] = 0;
//					ionTemperature[j] = Ti[0];
//					arealDensity[j] = ρR[0];
//				} else if (timeAxis[j] > time[time.length-1]) {
////					neutronYield[j] = 0;
//					ionTemperature[j] = Ti[time.length-1];
//					arealDensity[j] = ρR[time.length-1];
//				}
//				else {
////					neutronYield[j] = NumericalMethods.interp(timeAxis[j], time, Yn);
//					ionTemperature[j] = NumericalMethods.interp(timeAxis[j], time, Ti);
//					arealDensity[j] = NumericalMethods.interp(timeAxis[j], time, ρR);
//				}
//			}
//			double truePosterior = this.logPosterior(spectrum, neutronYield, ionTemperature, electronTemperature, bulkFlowVelocity, arealDensity, left, rite, false, finalSmoothing);
//			double trueBayesor = this.logPosterior(spectrum, neutronYield, ionTemperature, electronTemperature, bulkFlowVelocity, arealDensity, left, rite, false, 0);
//			System.out.printf("it got to %g, but the true anser is actually at %g.  in terms of likelihood alone, that's %g and %g.\n", ultimatePosterior, truePosterior, ultimateBayesor, trueBayesor);
//		}
//		catch (IOException e) {
//			System.err.println("the test thing didn't work.");
//		}

//		System.out.println("it arrived at");
//		System.out.println(Arrays.toString(ionTemperature));
//		System.out.println()
//		System.out.println("but wouldn't ")

		this.neutronYield = new Quantity[M];
		this.ionTemperature = new Quantity[M];
		this.arealDensity = new Quantity[M];
		for (int j = 0; j < timeAxis.length; j ++) {
			this.neutronYield[j] = new Quantity(
				  neutronYield[j], j, N*M);
			this.ionTemperature[j] = new Quantity(
				  ionTemperature[j], M+j, N*M);
			this.arealDensity[j] = new Quantity(
				  arealDensity[j], N*M);//, 2*M+j, N*M);
		}

		this.fitNeutronSpectrum = SpectrumGenerator.generateSpectrum( // and then interpret it
				neutronYield, ionTemperature, electronTemperature,
				bulkFlowVelocity, arealDensity, energyBins, timeBins);
		this.fitDeuteronSpectrum = this.response(
			  energyBins, timeBins, fitNeutronSpectrum, false, false);
		
		Quantity iBT = NumericalMethods.quadargmax(this.neutronYield); // index of max yield
		Quantity bangTime = NumericalMethods.interp(timeAxis, iBT); // time of max yield

		if (errorMode == ErrorMode.HESSIAN) { // calculate the error bars using second-derivatives
			logger.log(Level.INFO, "calculating error bars.");
			final int m = rite - left;
			double[] x0 = new double[N*m], dx = new double[N*m];
			for (int j = left; j < rite; j ++) { // set up for a finite-difference hessian calculation
				x0[      j - left] = neutronYield[j];
				dx[      j - left] = NumericalMethods.max(neutronYield)*1e-4;
				x0[  m + j - left] = ionTemperature[j];
				dx[  m + j - left] = 1e-2;
				x0[2*m + j - left] = arealDensity[j];
				dx[2*m + j - left] = 1e-3;
			}
			final double[] testNeutronYield = Arrays.copyOf(neutronYield, M);
			final double[] testIonTemperature = Arrays.copyOf(ionTemperature, M);
			final double[] testArealDensity = Arrays.copyOf(arealDensity, M);
			double[][] hessian = NumericalMethods.hessian((double[] x) -> {
				assert x.length == N*m;
				System.arraycopy(x, 0  , testNeutronYield, left, m);
				System.arraycopy(x, 1*m, testIonTemperature, left, m);
				System.arraycopy(x, 2*m, testArealDensity, left, m);
				return this.logPosterior(spectrum, testNeutronYield,
										 testIonTemperature, electronTemperature,
										 bulkFlowVelocity, testArealDensity,
										 active, false, finalSmoothing);
			}, x0, dx); // the do that
			NumericalMethods.coercePositiveSemidefinite(hessian);
			double[][] activeCovarianceMatrix = NumericalMethods.pseudoinv(hessian); // invert it to get the covariance
			for (int i = 0; i < hessian.length; i ++) {
				if (activeCovarianceMatrix[i][i] < 1/hessian[i][i]) // these are all approximations, and sometimes they violate the properties of positive semidefiniteness
					activeCovarianceMatrix[i][i] = 1/hessian[i][i]; // do what you must to make it work
			}
			covarianceMatrix = new double[N*M][N*M]; // then expand it to include the inactive dimensions
			for (int k1 = 0; k1 < N; k1 ++)
				for (int j1 = left; j1 < rite; j1++)
					for (int k2 = 0; k2 < N; k2++)
						for (int j2 = left; j2 < rite; j2++)
							covarianceMatrix[k1*M + j1][k2*M + j2] = activeCovarianceMatrix[k1*m + j1 - left][k2*m + j2 - left];
//			NumericalMethods.coerceSymmetric(covarianceMatrix);
			NumericalMethods.coercePositiveSemidefinite(covarianceMatrix);
//			System.out.println(Arrays.deepToString(covarianceMatrix));
		}
		else if (errorMode == ErrorMode.STATISTICS) {
			covarianceMatrix = new double[N*M][N*M];
			for (int j = 0; j < rite - left; j ++) {
				double primaryStatistics = 1 + neutronYield[j]*timeStep*1e15*this.gain(); // total deuteron yield from this neutron time bin
				double dsStatistics = 1 + neutronYield[j]*timeStep*1e15*this.gain()*arealDensity[j]/21.;
				covarianceMatrix[    j][    j] = Math.pow(neutronYield[j], 2)/primaryStatistics;
				covarianceMatrix[  M+j][  M+j] = Math.pow(ionTemperature[j], 2)*2/(primaryStatistics - 1);
				covarianceMatrix[2*M+j][2*M+j] = Math.pow(arealDensity[j], 2)/dsStatistics;
			}
		}
		else {
			covarianceMatrix = new double[N*M][N*M];
		}

		long endTime = System.currentTimeMillis();
		logger.info(String.format(Locale.US, "completed in %.2f minutes.",
				(endTime - startTime)/60000.));
		
		Quantity[] V = new Quantity[M]; // this is not really volume; it is (ρR)^(-2/3)
		for (int j = 0; j < V.length; j ++)
			V[j] = this.arealDensity[j].pow(-1.5);
		
		Quantity iPC = NumericalMethods.quadargmax(left, rite, this.arealDensity); // index of max compression
		Quantity peakCompress = NumericalMethods.interp(timeAxis, iPC); // time of max compression
		Quantity iTPeak = NumericalMethods.quadargmax(left, rite, this.ionTemperature);
		Quantity[] moments = new Quantity[5];
		for (int k = 0; k < moments.length; k ++)
			moments[k] = NumericalMethods.moment(k, timeBins, this.neutronYield);
		
		Quantity[] res = {
			  new Quantity((endTime - startTime)/1000., covarianceMatrix.length),
			  moments[0].times(timeStep), bangTime,
			  moments[2].sqrt().times(2.355), moments[3], moments[4],
			  peakCompress.minus(bangTime),
			  NumericalMethods.average(this.ionTemperature, this.neutronYield, left, rite),
			  NumericalMethods.quadInterp(this.ionTemperature, iTPeak),
			  NumericalMethods.quadInterp(this.ionTemperature, iBT),
			  NumericalMethods.derivative(timeAxis, this.ionTemperature, bangTime, .12, 1),
			  NumericalMethods.derivative(timeAxis, this.ionTemperature, bangTime, .12, 2),
			  NumericalMethods.average(this.arealDensity, this.neutronYield, left, rite),
			  NumericalMethods.quadInterp(this.arealDensity, iPC),
			  NumericalMethods.quadInterp(this.arealDensity, iBT),
			  NumericalMethods.derivative(timeAxis, this.arealDensity, bangTime, .12, 1),
			  NumericalMethods.derivative(timeAxis, V, bangTime, .12, 2).over(NumericalMethods.quadInterp(V, iBT)),
		}; // collect the figures of merit
		
		logger.info(String.format("Total yield (μ0):  %s", res[1].times(1e15).toString(covarianceMatrix)));
		logger.info(String.format("Bang time:         %s ns", res[2].toString(covarianceMatrix)));
		logger.info(String.format("Burn width (μ2):   %s ps", res[3].over(1e-3).toString(covarianceMatrix)));
		logger.info(String.format("Burn skewness (μ3):%s", res[4].toString(covarianceMatrix)));
		logger.info(String.format("Burn kurtosis (μ4):%s", res[5].toString(covarianceMatrix)));
		logger.info(String.format("Peak compression:  %s ps + BT", res[6].over(1e-3).toString(covarianceMatrix)));
		logger.info(String.format("Burn-averaged Ti:  %s keV", res[7].toString(covarianceMatrix)));
		logger.info(String.format("Ti at peak:        %s keV", res[8].toString(covarianceMatrix)));
		logger.info(String.format("Ti at BT:          %s keV", res[9].toString(covarianceMatrix)));
		logger.info(String.format("dTi/dt at BT:      %s keV/(100 ps)", res[10].over(1e1).toString(covarianceMatrix)));
		logger.info(String.format("d^2Ti/dt^2 at BT:  %s keV/ns^2", res[11].toString(covarianceMatrix)));
		logger.info(String.format("Burn-averaged \u03C1R:  %s g/cm^2", res[12].toString(covarianceMatrix)));
		logger.info(String.format("\u03C1R at peak:        %s g/cm^2", res[13].toString(covarianceMatrix)));
		logger.info(String.format("\u03C1R at BT:          %s g/cm^2", res[14].toString(covarianceMatrix)));
		logger.info(String.format("d\u03C1R/dt at BT:      %s g/cm^2/(100 ps)", res[15].over(1e1).toString(covarianceMatrix)));
		logger.info(String.format("d^2V/dt^2/V at BT: %s 1/ns^2", res[16].toString(covarianceMatrix)));
		return res;
	}

	/**
	 * a utility function for the fitting: get the error in this spectrum fit
	 * @param active whether to apply the prior to each timestep
	 * @param sumInEnergy whether to combine energies when comparing
	 * @return the log of the inverse posterior value for this spectrum
	 */
	private double logPosterior(double[][] spectrum,
								double[] neutronYield,
								double[] ionTemperature,
								double[] electronTemperature,
								double[] bulkFlowVelocity,
								double[] arealDensity,
								boolean[] active,
								boolean sumInEnergy,
								double smoothing) {
		for (int j = 0; j < timeAxis.length; j ++)
			if (Double.isNaN(neutronYield[j]))
				throw new IllegalArgumentException("you shouldn't be passing nan; whence did it come");
		for (int j = 0; j < timeAxis.length; j ++)
			if (neutronYield[j] < 0 || ionTemperature[j] < 0 ||
				  electronTemperature[j] < 0 || bulkFlowVelocity[j] < 0 ||
				  arealDensity[j] < 0)
				return Double.POSITIVE_INFINITY;

		double[][] neutrons = SpectrumGenerator.generateSpectrum(
			  neutronYield, ionTemperature, electronTemperature, bulkFlowVelocity, arealDensity,
			  energyBins, timeBins); // generate the neutron spectrum based on those
		double[][] signal = this.response(
			  energyBins, timeBins, neutrons, false, false);

		double gain = detector.gain();
		double totalError = 0; // negative log likelihood
		for (int j = 0; j < timeAxis.length; j ++) {
			int numEnergies = (sumInEnergy) ? 1 : spectrum.length;
			double[] experValues = new double[numEnergies];
			double[] theorValues = new double[numEnergies];
			double[] variances = new double[numEnergies];
			double[] backgrounds = new double[numEnergies];
			for (int i = 0; i < spectrum.length; i ++) {
				double energy = (energyBins[i] + energyBins[i+1])/2.;
				int index = (sumInEnergy) ? 0 : i;
				experValues[index] += spectrum[i][j];
				theorValues[index] += signal[i][j];
				variances[index] += detector.noise(energy, energyBins, timeBins);
				backgrounds[index] = detector.background(energy, energyBins, timeBins);
			}
			for (int i = 0; i < numEnergies; i ++) {
				double theorNumber = (theorValues[i] - backgrounds[i])/gain;
				double experNumber = (experValues[i] - backgrounds[i])/gain;
				if (variances[i] > 0) { // if this detector has significant noise
					double variance = variances[i] + theorNumber*gain*gain; // include pre-amplification poisson noise
					totalError += Math.pow(experValues[i] - theorValues[i], 2)/
						  (2*variance); // and use a Gaussian approximation
				}
				else { // if the detector noise is zero
					if (theorNumber > 0) {
						totalError += theorNumber - experNumber*Math.log(theorNumber); // use the exact Poisson distribution
					}
					else if (theorNumber == 0) {
						if (experNumber == 0)
							totalError += 0;
						else
							return Double.POSITIVE_INFINITY;
					}
					else {
						throw new IllegalArgumentException("What to do when expected "+theorNumber+" and observed "+experNumber+"?");
					}
				}
			}
		}

		double totalPenalty = 0; // negative log prior
		int bangIndex = NumericalMethods.argmax(neutronYield); // approximate BT locacion
		for (int j = 1; j < timeAxis.length; j ++) {
			double Y = (neutronYield[j] + neutronYield[j-1])/2;
			double Yp = (neutronYield[j] - neutronYield[j-1])/timeStep;
			if (j <= bangIndex) Yp *= -1;
			double z = (Y > 0) ? Yp/Y*1. : 0;
			if (z < 0) totalPenalty += .05*Math.exp(z);
			else       totalPenalty += .05*(1 + z + z*z/2.); // encourage a monotonically increasing yield before BT
		}
		for (int j = 0; j < timeAxis.length; j ++)
			totalPenalty += arealDensity[j]/5.0; // exponential prior on areal density
		for (int j = 0; j < timeAxis.length; j ++)
			totalPenalty += Math.pow((Math.log(ionTemperature[j]) - 0.2)/2.5, 2);
		for (double[] x: new double[][] {ionTemperature, arealDensity}) {
			for (int j = 0; j < timeAxis.length - 3; j ++) {
				if (active == null || (active[j] && active[j+1] && active[j+2] && active[j+3])) {
					double Ψpp = (x[j] - 3*x[j+1] + 3*x[j+2] - x[j+3])/
						  Math.pow(timeStep, 3);
					double Ψ = (x[j] + x[j+1] + x[j+2] + x[j+3])/4;
					if (Ψpp != 0)
						totalPenalty += smoothing*1e-10*Math.pow(Ψpp/Ψ, 2); // encourage a smooth Ti and ρR
				}
			}
		}

		return totalError + totalPenalty;
	}
	
	/**
	 * optimize certain elements of the input vector to maximize this function output.
	 * this routine is bizarrely convoluted, but if I don't do all of this it doesn't converge completely, and I don't
	 * understand why.
	 * @param func the objective function to optimize
	 * @param totalGuess the initial guess
	 * @param totalScale the scale lengths of the variables
	 * @param totalLower the lower bounds of the variables (only sometimes enforced)
	 * @param totalUpper the upper bounds of the variables (only sometimes enforced)
	 * @param activeDimensions which components of the state vector are allowed to be changed
	 */
	private double[] optimize(Function<double[], Double> func, double[] totalGuess,
							  double[] totalScale, double[] totalLower, double[] totalUpper,
							  boolean... activeDimensions) {
		if (totalGuess.length != totalScale.length)
			throw new IllegalArgumentException("Scale and guess must have the same length");
		
		final boolean[] active = new boolean[totalScale.length];
		for (int i = 0; i < active.length; i ++) { // expand the selected dimensions into a full array
			if (activeDimensions.length == 0)
				active[i] = true;
			else
				active[i] = activeDimensions[i%activeDimensions.length];
		}
		
		int numTotal = active.length; // count them
		int numActive = 0;
		for (boolean activeDimension: active)
			if (activeDimension)
				numActive ++;
		
		double[] activeGuess = new double[totalGuess.length/numTotal*numActive]; // so that you can make these arrays
		double[] activeScale = new double[totalScale.length/numTotal*numActive];
		double[] activeLower = new double[totalScale.length/numTotal*numActive];
		double[] activeUpper = new double[totalScale.length/numTotal*numActive];
		{ // this extra scope is here so I can redeclare j later
			int j = 0;
			for (int i = 0; i < totalGuess.length; i ++) {
				if (totalGuess[i] < totalLower[i])
					throw new IllegalArgumentException("initial guess ("+totalGuess[i]+") below bound at index "+i+" ("+totalLower[i]+", "+totalUpper[i]+")");
				if (totalGuess[i] > totalUpper[i])
					throw new IllegalArgumentException("initial guess ("+totalGuess[i]+") above bound at index "+i+" ("+totalLower[i]+", "+totalUpper[i]+")");
				if (active[i%active.length]) {
					activeGuess[j] = totalGuess[i];
					activeScale[j] = totalScale[i];
					activeLower[j] = totalLower[i];
					activeUpper[j] = totalUpper[i];
					j ++;
				}
			}
		}
		
		double[] totalParams = totalGuess.clone();
		logger.log(Level.FINER, Double.toString(func.apply(totalGuess)));
		activeGuess = Optimization.minimizeLBFGSB(
			  (activeParams) -> { // when you optimize
					updateArray(totalParams, activeParams, active); // only optimize a subset of the dimensions
					return func.apply(totalParams);
				},
			  activeGuess,
			  activeScale,
			  activeLower,
			  activeUpper,
			  0, 1.0*1e-2);
		updateArray(totalParams, activeGuess, active);
		
		double oldPosterior = Double.POSITIVE_INFINITY, newPosterior = func.apply(totalParams);
		logger.log(Level.FINER, Double.toString(newPosterior));
		MultivariateOptimizer optimizer = new PowellOptimizer(1e-14, 0.1);
		while (oldPosterior - newPosterior > 0.1) { // optimize it over and over; you'll get there eventually
			activeGuess = optimizer.optimize(
					GoalType.MINIMIZE,
					new ObjectiveFunction((activeParams) -> { // when you optimize
						updateArray(totalParams, activeParams, active);
						return func.apply(totalParams);
					}),
					new InitialGuess(activeGuess),
					new MultiDirectionalSimplex(activeScale),
					new MaxIter(10000),
					new MaxEval(100000)).getPoint();
			
			oldPosterior = newPosterior;
			updateArray(totalParams, activeGuess, active);
			newPosterior = func.apply(totalParams);
			logger.log(Level.FINER, Double.toString(newPosterior));
		}
		return totalParams;
	}
	
	private void updateArray(double[] totalVals, double[] activeVals, boolean[] active) {
		int j = 0;
		for (int i = 0; i < totalVals.length; i ++) {
			if (active[i%active.length]) {
				totalVals[i] = activeVals[j]; // you're only changing a subset of the parameters
				j ++;
			}
		}
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
	
	/**
	 * get the time bin centers in [ns]
	 */
	public double[] getTimeAxis() {
		return this.timeAxis;
	}
	
	/**
	 * get the energy bin centers in [keV]
	 */
	public double[] getEnergyAxis() {
		return this.energyAxis;
	}
	
	/**
	 * get the neutron yield mean values in [10^15/ns]
	 */
	public double[] getNeutronYield() {
		return NumericalMethods.modes(this.neutronYield);
	}
	
	/**
	 * get the neutron yield error bars in [10^15/ns]
	 */
	public double[] getNeutronYieldError() {
		return NumericalMethods.stds(this.neutronYield, this.covarianceMatrix);
	}
	
	/**
	 * get the ion temperature mean values in [keV]
	 */
	public double[] getIonTemperature() {
		return NumericalMethods.modes(this.ionTemperature);
	}
	
	/**
	 * get the ion temperature error bars in [keV]
	 */
	public double[] getIonTemperatureError() {
		return NumericalMethods.stds(this.ionTemperature, this.covarianceMatrix);
	}

	/**
	 * get the ρR mean values in [g/cm^2]
	 */
	public double[] getArealDensity() {
		return NumericalMethods.modes(this.arealDensity);
	}
	
	/**
	 * get the ρR error bars in [g/cm^2]
	 */
	public double[] getArealDensityError() {
		return NumericalMethods.stds(this.arealDensity, this.covarianceMatrix);
	}
	
	
	private static String[] appendErrorsToHeader() {
		String[] out = new String[2*HEADERS.length + 4];
		for (int i = 0; i < HEADERS.length; i ++) {
			if (i < 4)
				out[i] = HEADERS[i];
			else {
				out[2*(i-4)+4] = HEADERS[i];
				out[2*(i-4)+5] = HEADERS[i] + " error";
			}
		}
		return out;
	}
	
	
	/**
	 * different methods for computing error bars: assume they are all zero, estimate a
	 * hessian using finite differences and invert it into a covariance matrix, or calculate
	 * them directly with statistics.
	 * @author Justin Kunimune
	 *
	 */
	public enum ErrorMode {
		NONE, HESSIAN, STATISTICS
	}
	
	
//	public static final void main(String[] args) {
//		int n = 72;
//		double[] E = new double[n+1];
//		for (int i = 0; i <= n; i ++)
//			E[i] = 12 + i*4./n;
//		double[] Y = generateSpectrum(100, 8, 8, 0, 1, 0, E, false);
//		
//		double ds = 0;
//		for (int i = 0; i < n; i ++)
//			if (E[i] < 13)
//				ds += Y[i];
//		
//		System.out.println(Arrays.toString(E));
//		System.out.println(Arrays.toString(Y));
//		System.out.println(NumericalMethods.sum(Y));
//		System.out.println(ds/NumericalMethods.sum(Y));
//	}
}
