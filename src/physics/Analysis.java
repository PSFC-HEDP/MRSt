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
		"Burn-average \u03C1R (g/cm^2)", "\u03C1R at stagnation (g/cm^2)",
		"\u03C1R at BT (g/cm^2)", "d\u03C1R/dt at BT (g/cm^2/ns)", "d^2V/dt^2/V at BT (1/ns^2)",
		"Burn-average Ti (keV)", "Peak Ti (keV)",
		"Ti at BT (keV)", "dTi/dt at BT (keV/ns)", "d^2Ti/dt^2 at BT (keV/ns^2)",
		"Burn-average vi (km/s)", "dvi/dt at BT (km/s/ns)",
		}; // the names, units, and order of time-dependent burn parameters
	public static final String[] HEADERS_WITH_ERRORS = appendErrorsToHeader();

	public static final Random RANDOM = new Random(0);

	private static final double MIN_E = 12, MAX_E = 16; // histogram bounds [MeV]
	private static final int BUFFER = 4; // empty pixels to include simulate on each side [ns]
	private static final double E_BIN = .09, T_BIN = 30e-3; // bin sizes [MeV], [ns]

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
	private Quantity[][] measurements; // yield, ion temperature, electron temperature, velocity, areal density, and P2 asymmetry
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
			COSYMapping cosyMapping, double focalTilt,
			double precision, Logger logger) throws IOException {

		this.ionOptics = new IonOptics(
				foilDistance, foilWidth, foilHeight, foilThickness,
				apertureDistance, apertureWidth, apertureHeight,
				MIN_E, MAX_E, cosyMapping, focalTilt);
//		this.detector = new PulseDilationDriftTube(
//				ion, SUBSTRATE_THICKNESS, PHOTOCATHODE_THICKNESS,
//				focalTilt, PDDT_BIAS, MESH_LENGTH, DRIFT_LENGTH,
//				TIME_DILATION, MCT_POROSITY, MCT_GAIN, 100);
//		this.detector = new StreakCameraArray(
//			  2.5e-2, 200e-6, 5e-2,
//			  cosyMapping, apertureDistance, apertureHeight);
		this.detector = new PerfectDetector();

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
	 * @param errorBars should error bars be computed (it's rather intensive)?
	 * @return {computation time, Yn, BT, BW, skewness, kurtosis, peak compression, rho R (BT),
	 * \< rho R \>, d(rho R)/dt, \<Ti\>, dTi/dt, \<vi\>, dvi/dt} or null if it can't even
	 */
	private Quantity[] analyze(double[][] spectrum, ErrorMode errorBars) {
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

//		double[][] D = new double[energyBins.length-1][timeBins.length-1];
//		for (int i = 0; i < energyBins.length-1; i ++)
//			for (int j = 0; j < timeBins.length-1; j ++)
//				D[i][j] = Math.max(1, spectrum[i][j]);

//		double[][] gelf = Optimization.optimizeGelfgat(spectrum, D, this.rongTransferMatrix, 1e5);
		double[][] gelf = new double[spectrum.length][spectrum[0].length]; // start with the most basick inicial gess
		for (int i = 0; i < spectrum.length; i ++)
			for (int j = 0; j < spectrum[i].length; j ++)
				gelf[i][j] = spectrum[i][j]/this.gain(); // the spectrum scaled by the efficiency of the system
		
		double[] yieldGess = new double[timeAxis.length];
		for (int j = 0; j < timeAxis.length; j ++) {
			yieldGess[j] = 0;
			for (int i = 3; i < energyBins.length-1; i ++) // I'm not sure why the bottom few rows are so unusable
				yieldGess[j] += gelf[i][j];
		}
		
		final int bangIndex = NumericalMethods.argmax(yieldGess); // approximate BT locacion
		final int left = NumericalMethods.firstLocalMin(yieldGess);
		final int rite = NumericalMethods.lastLocalMin(yieldGess) + 1;
		
		double[] opt = new double[4*(rite - left)];
		double[] dimensionScale = new double[opt.length];
		double[] lowerBound = new double[opt.length];
		double[] upperBound = new double[opt.length];
		for (int j = 0; j < rite - left; j ++) {
			yieldGess[left+j] = Math.max(yieldGess[left+j], 1/this.gain());
			double scaledYieldGess = yieldGess[left+j]/timeStep/1e15;
			double[] gesses = {          scaledYieldGess,  4,    0, 1 }; //inicial gess
			double[] scales = {      0.5*scaledYieldGess, 10,  200, 1 }; // rough ranges of these variables
			double[] lowers = {                       0.,  0, -250, 0 }; // lower bounds
			double[] uppers = { Double.POSITIVE_INFINITY, 20,  250, 4 }; // upper bounds
			for (int k = 0; k < scales.length; k ++) {
				opt[4*j+k] = gesses[k];
				dimensionScale[4*j+k] = scales[k];
				lowerBound[4*j+k] = lowers[k];
				upperBound[4*j+k] = uppers[k];
			}
		}
		
		double[] smoothingRapper = new double[1];
		double[] baselineRapper = new double[1];
		
		double[] electronTemperature = new double[timeAxis.length];
		for (int j = 0; j < timeAxis.length; j ++)
			electronTemperature[j] = 3;
		
		Function<double[], Double> logPosterior = (double[] x) -> {
			final double smoothing = smoothingRapper[0];
			final double baseline = baselineRapper[0];
			
			double[][] params = new double[4][timeAxis.length];
			for (int j = 0; j < rite - left; j ++) {
				for (int k = 0; k < params.length; k ++) {
					params[k][left+j] = x[4*j+k]; // first unpack the state vector
					
					if (params[k][left+j] < lowerBound[4*j+k] ||
							params[k][left+j] > upperBound[4*j+k]) // check for illegal (prior = 0) values
						return Double.POSITIVE_INFINITY;
				}
			}
			
			double[][] teoSpectrum = SpectrumGenerator.generateSpectrum(
					params[0], params[1], electronTemperature, params[2], params[3],
					energyBins, timeBins); // generate the neutron spectrum based on those
			double[][] fitSpectrum = this.response(
				  energyBins, timeBins, teoSpectrum, false, false); // blur it according to the transfer matrix
			
			double error = 0; // negative Bayes factor (in nepers)
			for (int i = 0; i < spectrum.length; i ++) {
				for (int j = 0; j < spectrum[i].length; j ++) { // compute the error between it and the actual spectrum
					if (fitSpectrum[i][j] == 0) {
						if (spectrum[i][j] == 0)
							error += 0;
						else
							return Double.POSITIVE_INFINITY;
					}
					else if (fitSpectrum[i][j] > 0) {
						error += fitSpectrum[i][j] - spectrum[i][j]*Math.log(fitSpectrum[i][j]);
//						double θ = fitSpectrum[i][j];
//						double ηe = detector.efficiency(energyAxis[i]);
//						double log_p0 = θ*(Math.exp(-ηe) - 1);
//						if (spectrum[i][j] == 0) {
//							error += -log_p0;
//						}
//						else if (spectrum[i][j] > 0) {
//							double n = spectrum[i][j];
//							double K = detector.gain();
//							double μ = K*ηe*θ;
//							if (ηe*θ > 100) {
//								μ = θ/(1 - Math.exp(-θ));
//								μ = K*ηe*μ/(1 - Math.exp(-ηe*μ));
//							}
//							double β = Math.max(1/μ, (.32/Math.sqrt(ηe) + (.09+.01*ηe)/Math.pow(θ,.45))/K);
//							double α = β*μ;
//							double log_Γα = (α - 1 > 1e-50) ? (α - 1)*Math.log(α - 1) - (α - 1) : 0;
//							error += -Math.log(1 - Math.exp(log_p0)) - (α - 1)*Math.log(n) + β*n - α*Math.log(β) + log_Γα;
//						}
//						else
//							throw new IllegalArgumentException("What to do when expected "+fitSpectrum[i][j]+" and observed "+spectrum[i][j]+"?");
					}
				}
			}
			
			double penalty = baseline; // negative log of prior (ignoring global normalization)
			for (int j = left; j < rite; j ++) {
				penalty += Math.pow(params[2][j]/50, 2)/2; // gaussian prior on velocity
//				penalty += params[1][j]/10.0 - Math.log(params[1][j])/2.; // gamma prior on temp
				penalty += params[3][j]/2.0; // exponential prior on areal density
			}
			
			for (int j = left + 1; j < rite; j ++) {
				double Y = (params[0][j] + params[0][j-1])/2;
				double Yp = (params[0][j] - params[0][j-1])/timeStep;
				if (j <= bangIndex) Yp *= -1;
				if (Y > 0) {
					double z = Yp/Y*1.;
					if (z < 0) penalty += smoothing/20*Math.exp(z);
					else       penalty += smoothing/20*(1 + z + z*z/2.); // encourage a monotonically increasing yield before BT
				}
			}
			
			for (int k = 1; k <= 3; k += 2) {
				for (int j = left; j < rite - 3; j ++) {
					double Ψpp = (params[k][j] - 3*params[k][j+1] + 3*params[k][j+2] - params[k][j+3])/
							Math.pow(timeStep, 3);
					double Ψ = (params[k][j] + params[k][j+1] + params[k][j+2] + params[k][j+3])/4;
					if (Ψpp != 0)
						penalty += smoothing*5e-10*Math.pow(Ψpp/Ψ, 2); // encourage a smooth Ti and ρR
				}
			}
			
			for (int j = left + 1; j < rite - 1; j ++) {
				double Vpp = (params[2][j-1] - 2*params[2][j] + params[2][j+1])/
						Math.pow(timeStep, 3);
				penalty += smoothing*2e-6*Math.pow(Vpp/50, 2); // encourage a smooth ion velocity
			}
			
//			System.out.println(penalty+" + "+error);
			if (baseline == 0) { // the first time this function is called
				baselineRapper[0] = -(penalty + error); // save its value as the baseline
				return 0.0;
			}
			return penalty + error;
		};
		
		smoothingRapper[0] = 100;
		logger.log(Level.FINE, "Performing ruff fit pass...");
		opt = optimize(logPosterior, opt, dimensionScale, lowerBound, upperBound, 10.*precision);
		smoothingRapper[0] = 10;
		logger.log(Level.FINE, "Performing medium fit pass...");
		opt = optimize(logPosterior, opt, dimensionScale, lowerBound, upperBound, 1.*precision);
		smoothingRapper[0] = 3;
		logger.log(Level.FINE, "Performing careful fit pass...");
		opt = optimize(logPosterior, opt, dimensionScale, lowerBound, upperBound, .1*precision);
		smoothingRapper[0] = 1;
		logger.log(Level.FINE, "Performing final fit pass...");
		opt = optimize(logPosterior, opt, dimensionScale, lowerBound, upperBound, .001*precision);

		this.measurements = new Quantity[4][timeAxis.length]; // unpack the optimized vector
		for (int j = 0; j < timeAxis.length; j ++) {
			if (j >= left && j < rite) {
				for (int k = 0; k < measurements.length; k ++) {
					double[] grad = new double[opt.length];
					grad[4*(j-left)+k] = 1;
					measurements[k][j] = new Quantity(opt[4*(j-left)+k], grad);
				}
			}
			else {
				measurements[0][j] = new Quantity(0, opt.length);
				for (int k = 1; k < measurements.length; k ++)
					measurements[k][j] = new Quantity(Double.NaN, opt.length);
			}
		}
		
		this.fitNeutronSpectrum = SpectrumGenerator.generateSpectrum( // and then interpret it
				getNeutronYield(), getIonTemperature(), electronTemperature,
				getFlowVelocity(), getArealDensity(), energyBins, timeBins);
		this.fitDeuteronSpectrum = this.response(
			  energyBins, timeBins, fitNeutronSpectrum, false, false);
		
		Quantity iBT = NumericalMethods.quadargmax(measurements[0]); // index of max yield
		Quantity bangTime = NumericalMethods.interp(timeAxis, iBT); // time of max yield
		
		if (errorBars == ErrorMode.HESSIAN) {
			logger.log(Level.INFO, "calculating error bars.");
			double c = logPosterior.apply(opt); // start by getting the actual value
			double[] step = new double[opt.length]; // and the values in all basis directions
			for (int i = 0; i < step.length; i ++) {
				double dxi = dimensionScale[i]*1e-4;
				opt[i] += dxi;
				step[i] = logPosterior.apply(opt);
				opt[i] -= dxi;
			}
			double[][] hessian = new double[opt.length][opt.length]; // then go for the second derivatives
			for (int i = 0; i < hessian.length; i ++) { // TODO: move this into a separate funccion
				double dxi = dimensionScale[i]*1e-4;
				double r = step[i];
				opt[i] -= dxi;
				double l = logPosterior.apply(opt);
				opt[i] += dxi;
				if (Double.isInfinite(l)) { // if we are at a bound
					hessian[i][i] = Math.pow((r - c)/dxi, 2); // approximate this exponential-ish distribution as gaussian
					for (int j = 0; j < i; j ++)
						hessian[i][j] = hessian[j][i] = 0; // and reset any diagonal terms that previously involved this
				}
				else {
					hessian[i][i] = (r - 2*c + l)/(dxi*dxi); // otherwise approximate it as gaussian
					for (int j = 0; j < opt.length; j ++) { // and get some diagonal terms
						double dxj = dimensionScale[j]*1e-4;
						double u = step[j];
						opt[j] += dxj;
						opt[i] += dxi;
						double ur = logPosterior.apply(opt);
						opt[i] -= dxi;
						opt[j] -= dxj;
						hessian[i][j] = hessian[j][i] = (ur - u - r + c)/(dxi*dxj);
					}
				}
			}
			for (int i = 0; i < hessian.length; i ++) {
				if (hessian[i][i] < 0)
					hessian[i][i] = Double.NaN;
			}
			NumericalMethods.coercePositiveSemidefinite(hessian);
			covarianceMatrix = NumericalMethods.pseudoinv(hessian);
			for (int i = 0; i < hessian.length; i ++) {
				if (covarianceMatrix[i][i] < 1/hessian[i][i]) // these are all approximations, and sometimes they violate the properties of positive semidefiniteness
					covarianceMatrix[i][i] = 1/hessian[i][i]; // do what you must to make it work
			}
//			NumericalMethods.coerceSymmetric(covarianceMatrix);
			NumericalMethods.coercePositiveSemidefinite(covarianceMatrix);
		}
		else if (errorBars == ErrorMode.STATISTICS) {
			covarianceMatrix = new double[opt.length][opt.length];
			for (int j = 0; j < rite - left; j ++) {
				double statistics = 1 + opt[4*j+0]*timeStep*1e15*this.gain(); // total deuteron yield from this neutron time bin
				double dsStatistics = 1 + opt[4*j+0]*timeStep*1e15*this.gain()*opt[4*j+3]/21.;
				covarianceMatrix[4*j+0][4*j+0] = Math.pow(opt[4*j+0], 2)/statistics;
				covarianceMatrix[4*j+1][4*j+1] = Math.pow(opt[4*j+1], 2)*2/(statistics - 1);
				covarianceMatrix[4*j+2][4*j+2] = .4034*14*opt[4*j+1]/1e3/statistics/Math.pow(.54e-3, 2);
				covarianceMatrix[4*j+3][4*j+3] = Math.pow(opt[4*j+3], 2)/dsStatistics;
			}
		}
		else {
			covarianceMatrix = new double[opt.length][opt.length];
		}
		
		long endTime = System.currentTimeMillis();
		logger.info(String.format(Locale.US, "completed in %.2f minutes.",
				(endTime - startTime)/60000.));
		
		Quantity[] V = new Quantity[measurements[3].length]; // this is not really volume; it is (ρR)^(-2/3)
		for (int j = 0; j < V.length; j ++)
			V[j] = measurements[3][j].pow(-1.5);
		
		Quantity iPC = NumericalMethods.quadargmax(left, rite, measurements[3]); // index of max compression
		Quantity peakCompress = NumericalMethods.interp(timeAxis, iPC); // time of max compression
		Quantity iTPeak = NumericalMethods.quadargmax(left, rite, measurements[1]);
		Quantity[] moments = new Quantity[5];
		for (int k = 0; k < moments.length; k ++)
			moments[k] = NumericalMethods.moment(k, timeBins, measurements[0]);
		
		Quantity[] res = {
				new Quantity((endTime - startTime)/1000., covarianceMatrix.length),
				moments[0].times(timeStep), bangTime,
				moments[2].sqrt().times(2.355), moments[3], moments[4],
				peakCompress.minus(bangTime),
				NumericalMethods.average(measurements[3], measurements[0], left, rite),
				NumericalMethods.quadInterp(measurements[3], iPC),
				NumericalMethods.quadInterp(measurements[3], iBT),
				NumericalMethods.derivative(timeAxis, measurements[3], bangTime, .12, 1),
				NumericalMethods.derivative(timeAxis, V, bangTime, .12, 2).over(NumericalMethods.quadInterp(V, iBT)),
				NumericalMethods.average(measurements[1], measurements[0], left, rite),
				NumericalMethods.quadInterp(measurements[1], iTPeak),
				NumericalMethods.quadInterp(measurements[1], iBT),
				NumericalMethods.derivative(timeAxis, measurements[1], bangTime, .12, 1),
				NumericalMethods.derivative(timeAxis, measurements[1], bangTime, .12, 2),
				NumericalMethods.average(measurements[2], measurements[0], left, rite),
				NumericalMethods.derivative(timeAxis, measurements[2], bangTime, .12, 1),
		}; // collect the figures of merit
		
		logger.info(String.format("Total yield (μ0):  %s", res[1].times(1e15).toString(covarianceMatrix)));
		logger.info(String.format("Bang time:         %s ns", res[2].toString(covarianceMatrix)));
		logger.info(String.format("Burn width (μ2):   %s ps", res[3].over(1e-3).toString(covarianceMatrix)));
		logger.info(String.format("Burn skewness (μ3):%s", res[4].toString(covarianceMatrix)));
		logger.info(String.format("Burn kurtosis (μ4):%s", res[5].toString(covarianceMatrix)));
		logger.info(String.format("Peak compression:  %s ps + BT", res[6].over(1e-3).toString(covarianceMatrix)));
		logger.info(String.format("Burn-averaged \u03C1R:  %s g/cm^2", res[7].toString(covarianceMatrix)));
		logger.info(String.format("\u03C1R at peak:        %s g/cm^2", res[8].toString(covarianceMatrix)));
		logger.info(String.format("\u03C1R at BT:          %s g/cm^2", res[9].toString(covarianceMatrix)));
		logger.info(String.format("d\u03C1R/dt at BT:      %s g/cm^2/(100 ps)", res[10].over(1e1).toString(covarianceMatrix)));
		logger.info(String.format("d^2V/dt^2/V at BT: %s 1/ns^2", res[11].toString(covarianceMatrix)));
		logger.info(String.format("Burn-averaged Ti:  %s keV", res[12].toString(covarianceMatrix)));
		logger.info(String.format("Ti at peak:        %s keV", res[13].toString(covarianceMatrix)));
		logger.info(String.format("Ti at BT:          %s keV", res[14].toString(covarianceMatrix)));
		logger.info(String.format("dTi/dt at BT:      %s keV/(100 ps)", res[15].over(1e1).toString(covarianceMatrix)));
		logger.info(String.format("d^2Ti/dt^2 at BT:  %s keV/ns^2", res[16].toString(covarianceMatrix)));
		logger.info(String.format("Burn-averaged vi:  %s km/s", res[17].toString(covarianceMatrix)));
		logger.info(String.format("dvi/dt at BT:      %s μm/ns/(100 ps)", res[18].over(1e1).toString(covarianceMatrix)));
		return res;
	}
	
	/**
	 * optimize certain elements of the input vector to maximize this function output.
	 * this routine is bizarrely convoluted, but if I don't do all of this it doesn't converge completely, and I don't
	 * understand why.
	 * @param func the objective function to optimize
	 * @param totalGuess the initial guess
	 * @param totalScale the scale lengths of the variables
	 * @param threshold the termination condition threshold
	 * @param activeDimensions which components of the state vector are allowed to be changed
	 */
	private double[] optimize(Function<double[], Double> func, double[] totalGuess,
			double[] totalScale, double[] totalLower, double[] totalUpper, double threshold,
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
					throw new IllegalArgumentException("initial guess below bound at index "+i);
				if (totalGuess[i] > totalUpper[i])
					throw new IllegalArgumentException("initial guess above bound at index "+i);
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
				0, threshold*1e-2);
		updateArray(totalParams, activeGuess, active);
		
		double oldPosterior = Double.POSITIVE_INFINITY, newPosterior = func.apply(totalParams);
		logger.log(Level.FINER, Double.toString(newPosterior));
		MultivariateOptimizer optimizer = new PowellOptimizer(1e-14, 1);
		while (oldPosterior - newPosterior > threshold) { // optimize it over and over; you'll get there eventually
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
		return NumericalMethods.modes(this.measurements[0]);
	}
	
	/**
	 * get the neutron yield error bars in [10^15/ns]
	 */
	public double[] getNeutronYieldError() {
		return NumericalMethods.stds(this.measurements[0], this.covarianceMatrix);
	}
	
	/**
	 * get the ion temperature mean values in [keV]
	 */
	public double[] getIonTemperature() {
		return NumericalMethods.modes(this.measurements[1]);
	}
	
	/**
	 * get the ion temperature error bars in [keV]
	 */
	public double[] getIonTemperatureError() {
		return NumericalMethods.stds(this.measurements[1], this.covarianceMatrix);
	}
	
	/**
	 * get the bulk flow mean values in [km/s]
	 */
	public double[] getFlowVelocity() {
		return NumericalMethods.modes(this.measurements[2]);
	}
	
	/**
	 * get the bulk flow mean values in [km/s]
	 */
	public double[] getFlowVelocityError() {
		return NumericalMethods.stds(this.measurements[2], this.covarianceMatrix);
	}
	
	/**
	 * get the ρR mean values in [g/cm^2]
	 */
	public double[] getArealDensity() {
		return NumericalMethods.modes(this.measurements[3]);
	}
	
	/**
	 * get the ρR error bars in [g/cm^2]
	 */
	public double[] getArealDensityError() {
		return NumericalMethods.stds(this.measurements[3], this.covarianceMatrix);
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
