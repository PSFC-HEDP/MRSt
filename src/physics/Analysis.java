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
import physics.Detector.DetectorConfiguration;
import physics.IonOptics.IonOpticConfiguration;
import util.COSYMapping;
import util.Math2;
import util.Math2.Quantity;
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
		  "Yield", "Computation time (s)",
		  "Total yield (10^15)", "Bang time (ns)", "Burn width (ns)",
		  "Burn skewness", "Burn kurtosis", "Stagnation - BT (ns)",
		  "Burn-average Ti (keV)", "Peak Ti (keV)",
		  "Ti at stagnation (keV)", "Ti at BT (keV)",
		  "dTi/dt at stagnation (g/cm^2/ns)", "dTi/dt at BT (keV/ns)", "d^2Ti/dt^2 at BT (keV/ns^2)",
		  "Burn-average \u03C1R (g/cm^2)", "\u03C1R at stagnation (g/cm^2)",
		  "\u03C1R at BT (g/cm^2)", "d\u03C1R/dt at BT (g/cm^2/ns)",
		  "d^2V/dt^2/V at BT (1/ns^2)",
		}; // the names, units, and order of time-dependent burn parameters
	public static final String[] HEADERS_WITH_ERRORS = appendErrorsToHeader();

	public static final Random MC_RANDOM = new Random(0);
	public static final Random NOISE_RANDOM = new Random(0);

	private static final double MIN_E = 12, MAX_E = 16; // histogram bounds [MeV]
	private static final int BUFFER = 5; // empty pixels to include simulate on each side [ns]

	public static final double DEFAULT_ENERGY_BIN = 50e-3;
	public static final double DEFAULT_TIME_BIN = 18e-3;
	private static final double DEFAULT_PRECISION = .01;

//	private static final double SUBSTRATE_THICKNESS = 100; // [μm]
//	private static final double PHOTOCATHODE_THICKNESS = .1; // [μm]
//	private static final double PDDT_BIAS = 1e3; // [V]
//	private static final double MESH_LENGTH = 1e-3; // [m]
//	private static final double DRIFT_LENGTH = 1e0; // [m]
//	private static final double TIME_DILATION = 20;
//	private static final double MCT_POROSITY = .70;
//	private static final double MCT_GAIN = 1e4;

	private static final double ELECTRON_TEMPERATURE = 4;
	private static final double BULK_FLOW_VELOCITY = 0;


	private final IonOptics ionOptics; // the ion optic system
	private final Detector detector; // the detector system

	private final double precision; // factor by which to ease the convergence conditions

	private final double[] energyBins; // endpoints of E bins for inferred spectrum (MeV n)
	private final double[] deuteronEnergyBins; // endpoints of E bins for inferred spectrum (MeV d)
	private double[] timeBins; // endpoints of time bins for inferred spectrum [ns]
	private double[][] deuteronSpectrum; // time-corrected deuteron counts
	private double[][] signalDistribution; // 2D-resolved signal measurement
	private double[][] idealSignalDistribution; // based on input neutron counts with no noise or gaps
	private double[][] fitNeutronSpectrum; // backward-fit neutron counts
	private double[][] fitDeuteronSpectrum; // backward-fit deuteron counts (this should be similar to deuteronSpectrum)
	private double[][] fitSignalDistribution; // backward-fit deuteron counts

	private final double[] energyAxis; // centers of energy bins (MeV neutron)
	private final double[] deuteronEnergyAxis; // centers of energy bins (MeV deuteron)
	private final double preferredTimeStep;
	private double timeStep;
	private double[] timeAxis; // centers of time bins (ns)
	private Quantity[] neutronYield; // 1e15/ns
	private Quantity[] ionTemperature; // keV
	private Quantity[] arealDensity; // g/cm^2
	private double[][] covarianceMatrix; // and covariances that go with all of these

	private final Logger logger; // for logging

	/**
	 * perform some preliminary calculations for the provided configuration.
	 * @throws IOException if one or more of the stopping power tables cannot be
	 *                     accessd for any reason
	 */
	public Analysis(
		  IonOpticConfiguration ionOpticConfiguration,
		  DetectorConfiguration detectorConfiguration,
		  double calibrationPrecision, boolean reuseMatrix,
		  Logger logger) throws IOException {
		this(new IonOptics(ionOpticConfiguration,
						   detectorConfiguration.cosyFile,
						   detectorConfiguration.tiltAngle,
						   detectorConfiguration.offset,
						   calibrationPrecision,
						   reuseMatrix),
			 detectorConfiguration,
			 logger);
	}

	/**
	 * perform some preliminary calculations for the provided configuration.
	 * @throws IOException if one or more of the stopping power tables cannot be
	 *                     accessd for any reason
	 */
	public Analysis(
		  IonOpticConfiguration ionOpticConfiguration,
		  DetectorConfiguration detectorConfiguration,
		  double calibrationPrecision, boolean reuseMatrix,
		  double eBin, double tBin, double analysisPrecision,
		  Logger logger) throws IOException {
		this(new IonOptics(ionOpticConfiguration,
						   detectorConfiguration.cosyFile,
						   detectorConfiguration.tiltAngle,
						   detectorConfiguration.offset,
						   calibrationPrecision, reuseMatrix),
			 detectorConfiguration,
			 eBin, tBin,
			 analysisPrecision, logger);
	}

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
	 * @throws IOException if one or more of the stopping power tables cannot be
	 *                     accessd for any reason
	 */
	public Analysis(
			double foilDistance, double foilWidth, double foilHeight, double foilThickness,
			double apertureDistance, double apertureWidth, double apertureHeight,
			COSYMapping cosyMapping,
			DetectorConfiguration detectorConfiguration,
			double calibrationPrecision, boolean reuseMatrix,
			Logger logger) throws IOException {

		this(new IonOptics(
				foilDistance, foilWidth, foilHeight, foilThickness,
				apertureDistance, apertureWidth, apertureHeight,
				MIN_E, MAX_E, cosyMapping, detectorConfiguration.tiltAngle,
				detectorConfiguration.offset, calibrationPrecision, reuseMatrix),
		     detectorConfiguration, logger);
	}

	/**
	 * perform some preliminary calculations for the provided configuration.
	 */
	public Analysis(
			IonOptics ionOptics, DetectorConfiguration detectorConfiguration,
			Logger logger) {
		this(ionOptics,
		     Detector.newDetector(detectorConfiguration, ionOptics),
		     logger);
	}

	/**
	 * perform some preliminary calculations for the provided configuration.
	 */
	public Analysis(
		  IonOptics ionOptics, DetectorConfiguration detectorConfiguration,
		  double eBin, double tBin, double analysisPrecision, Logger logger) {
		this(ionOptics,
		     Detector.newDetector(detectorConfiguration, ionOptics),
		     eBin, tBin, analysisPrecision, logger);
	}

	/**
	 * perform some preliminary calculations for the provided configuration.
	 */
	public Analysis(
			IonOptics ionOptics, Detector detector, Logger logger) {
		this(ionOptics, detector,
			 DEFAULT_ENERGY_BIN, DEFAULT_TIME_BIN, DEFAULT_PRECISION,
			 logger);
	}

	/**
	 * perform some preliminary calculations for the provided configuration.
	 */
	public Analysis(
		  IonOptics ionOptics, Detector detector,
		  double eBin, double tBin, double precision,
		  Logger logger) {

		this.ionOptics = ionOptics;
		this.detector = detector;

		this.precision = precision;

		this.energyBins = new double[(int) ((MAX_E - MIN_E)/eBin + 1)];
		this.deuteronEnergyBins = new double[energyBins.length];
		for (int i = 0; i < energyBins.length; i ++) {
			this.energyBins[i] = (MIN_E + i*(MAX_E - MIN_E)/(energyBins.length - 1));
			this.deuteronEnergyBins[i] = this.energyBins[i]*ionOptics.energyFactor;
		}

		this.energyAxis = new double[energyBins.length-1];
		this.deuteronEnergyAxis = new double[energyBins.length-1];
		for (int i = 0; i < energyBins.length-1; i ++) {
			this.energyAxis[i] = (this.energyBins[i] + this.energyBins[i + 1])/2;
			this.deuteronEnergyAxis[i] = this.energyAxis[i]*ionOptics.energyFactor;
		}

		this.preferredTimeStep = tBin;

		this.logger = logger;
	}

	/**
	 * establish the time bins and related quantities
	 * @param spectrumTimeBins the time bins of the input spectrum, to use as
	 *                         a reference
	 */
	private void instantiateTimeAxis(double[] spectrumTimeBins) {
		double timeOffset = NOISE_RANDOM.nextDouble()*preferredTimeStep;
		double minT = spectrumTimeBins[0] - BUFFER*preferredTimeStep + timeOffset;
		double maxT = spectrumTimeBins[spectrumTimeBins.length-1] + BUFFER*preferredTimeStep + timeOffset;
		this.timeBins = new double[(int) ((maxT - minT)/preferredTimeStep + 1)];
		for (int i = 0; i < timeBins.length; i ++)
			this.timeBins[i] = minT + i*(maxT - minT)/(timeBins.length-1);

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


	/**
	 * the detection efficiency of the complete system
	 * @param energy the energy at which to evaluate the efficiency
	 */
	public double efficiency(double energy) {
		return this.ionOptics.efficiency(energy)
				*this.detector.efficiency(energy, true);
	}


	/**
	 * the weighted average of the efficiency of the complete system
	 * @param temperature the ion temperature with which to weit it
	 * @param arealDensity the ρR with which to weit it
	 * @param minEnergy the lower neutron energy bound (MeV n)
	 * @param maxEnergy the upper neutron energy bound (MeV n)
	 */
	public double averageEfficiency(double temperature, double arealDensity,
									double minEnergy, double maxEnergy) {
		double[] weits = SpectrumGenerator.generateSpectrum(
				1, temperature, temperature, 0, arealDensity,
				energyBins, temperature == 0);
		double[] ionOpticEfficiency = ionOptics.response(energyBins, weits, false, false);
		double weitedEfficiency = 0;
		for (int i = 0; i < energyAxis.length; i ++) // perform a weited mean
			if (energyAxis[i] >= minEnergy && energyAxis[i] < maxEnergy)
				weitedEfficiency += ionOpticEfficiency[i]*detector.efficiency(energyAxis[i], true);
		return weitedEfficiency;
	}


	/**
	 * figure out what the nonstochastic response function is when you average
	 * over a range of energies (you know, like a low-res slit does)
	 * @param temperature the ion temperature with which to weit it
	 * @param arealDensity the ρR with which to weit it
	 * @param minEnergy the lower neutron energy bound (MeV n)
	 * @param maxEnergy the upper neutron energy bound (MeV n)
	 * @return a dimensionless kernel that you convolve with (n/bin) to get
	 * (signal/bin), without background
	 */
	public double[] averageResponse(double temperature, double arealDensity,
	                                double minEnergy, double maxEnergy) {
		double[] genericNeutronDistribution = SpectrumGenerator.generateSpectrum(
			  1, temperature, temperature, 0, arealDensity,
			  energyBins, temperature == 0);
		double[][] genericNeutronSpectrum = new double[energyAxis.length][timeAxis.length];
		for (int i = 0; i < energyAxis.length; i ++)
			genericNeutronSpectrum[i][timeAxis.length/2] = genericNeutronDistribution[i];
		double[][] typicalOpticsResponse = this.ionOptics.response(
				energyBins, timeBins, genericNeutronSpectrum,
				false, false);
		double[][] fullResponse = this.detector.response(
			  energyBins, timeBins, typicalOpticsResponse,
			  false, false, true);
		double[] bandResponse = new double[timeAxis.length];
		for (int i = 0; i < energyAxis.length; i ++)
			if (energyAxis[i] >= minEnergy && energyAxis[i] < maxEnergy)
				for (int j = 0; j < timeAxis.length; j ++)
					bandResponse[j] += fullResponse[i][j];
		return bandResponse;
	}


	/**
	 * compute the total amount of noise in this range
	 * @return the variance in this region (signal^2/timebin)
	 */
	public double totalVariance(double minEnergy, double maxEnergy) {
		double total = 0;
		for (double energy: energyAxis)
			if (energy >= minEnergy && energy <= maxEnergy)
				total += detector.noise(energy, energyBins, timeBins, true);
		return total;
	}


	/**
	 * compute the total amount of background in this range
	 * @return the background level in this region (signal/timebin)
	 */
	public double totalBackground(double minEnergy, double maxEnergy) {
		double total = 0;
		for (double energy: energyAxis)
			if (energy >= minEnergy && energy <= maxEnergy)
				total += detector.background(energy, energyBins, timeBins, true);
		return total;
	}


	/**
	 * compute the total response to an implosion with the given neutron spectrum
	 * and save and analyze it.
	 * @param energies the energies that describe the rows of counts [MeV]
	 * @param times the times that describe the columns of counts [ns]
	 * @param neutronSpectrum the number of neutrons in each time and energy bin.
	 *                        each row corresponds to one energy, and each column
	 *                        one time.
	 * @param errorBars whether to bother computing error bars
	 * @return {computation time, 0, Yn, err, BT, err, BW, err, skewness, err,
	 *          kurtosis, err, peak compression, err, rho R (BT), err,
	 *          \< rho R \>, err, d(rho R)/dt, err, \<Ti\>, err, dTi/dt, err,
	 *          \<vi\>, err, dvi/dt, err}
	 */
	public double[] respondAndAnalyze(double[] energies, double[] times,
	                                  double[][] neutronSpectrum,
	                                  ErrorMode errorBars) {
		if (neutronSpectrum.length != energies.length-1 || neutronSpectrum[0].length != times.length-1)
			throw new IllegalArgumentException("These dimensions don't make any sense.");

//		logger.info("responding to spectrum with yield "+
//				            Math2.sum(neutronSpectrum));

		this.instantiateTimeAxis(times); // first, create the detector time bins

		logger.info("beginning Monte Carlo computation.");
		long startTime = System.currentTimeMillis();

		// input neutron counts
		neutronSpectrum = Math2.downsample(
				times, energies, neutronSpectrum, timeBins, energyBins); // first rebin the spectrum
		this.deuteronSpectrum = this.ionOptics.response(
				energyBins, timeBins, neutronSpectrum,
				true, true);
		this.signalDistribution = this.detector.response(
			  energyBins, timeBins, deuteronSpectrum,
			  true, true, true);
		this.idealSignalDistribution = this.detector.response(
				energyBins, timeBins,
				this.ionOptics.response(
						energyBins, timeBins, neutronSpectrum,
						false, false),
				false, true, false);

		long endTime = System.currentTimeMillis();
		logger.info(String.format(Locale.US, "completed in %.2f minutes.",
		                          (endTime - startTime)/60000.));

		double saturation = Math2.max(signalDistribution)/
				this.detector.saturationLimit(14, energyBins, timeBins);
		if (saturation >= 1)
			logger.warning(String.format("The detector is saturating by about %.2f%%",
			                             saturation/1e-2));
		else
			logger.info(String.format("The detector is %.2f%% of the way to saturation",
									  saturation/1e-2));

		double totalEfficiency = this.averageEfficiency(7., .8, 12, 16);
		double opticsEfficiency = this.ionOptics.efficiency(14);
		System.out.printf("I think the ion-optic efficiency is %.3g and detector efficiency is %.3g\n", opticsEfficiency, totalEfficiency/opticsEfficiency);

		if (Math2.max(signalDistribution) == 0) {
			logger.log(Level.SEVERE, "There were no deuterons detected.");
			return null;
		}
		else {
			logger.log(Level.INFO, String.format("There were %.1g deuterons detected.", Math2.sum(signalDistribution)));
		}

		logger.info("beginning fit process.");
		startTime = System.currentTimeMillis();

		analyze(signalDistribution, errorBars);

		int dofs = covarianceMatrix.length;

		this.fitNeutronSpectrum = SpectrumGenerator.generateSpectrum(
				this.getNeutronYield(), this.getIonTemperature(),
				Math2.full(ELECTRON_TEMPERATURE, timeAxis.length),
				Math2.full(BULK_FLOW_VELOCITY, timeAxis.length),
				this.getArealDensity(),
				energyBins, timeBins);
		this.fitDeuteronSpectrum = this.ionOptics.response(
			  energyBins, timeBins, fitNeutronSpectrum,
			  false, false);
		this.fitSignalDistribution = this.detector.response(
			  energyBins, timeBins, fitDeuteronSpectrum,
			  false, true, false);

		endTime = System.currentTimeMillis();
		logger.info(String.format("completed in %.2f minutes.",
		                          (endTime - startTime)/60000.));

		Quantity[] V = new Quantity[timeAxis.length]; // this is not really volume; it is (ρR)^(-2/3)
		for (int j = 0; j < V.length; j ++)
			V[j] = this.arealDensity[j].pow(-1.5);

		Quantity iBT = Math2.quadargmax(this.neutronYield); // index of max yield
		Quantity bangTime = Math2.interp(timeAxis, iBT); // time of max yield

		Quantity peak = Math2.quadInterp(this.neutronYield, iBT);
		int left, rite;
		for (left = 0; left + 1 < iBT.value; left ++)
			if (neutronYield[left].value > peak.value*1e-2)
				break;
		for (rite = timeAxis.length; rite > iBT.value; rite --)
			if (neutronYield[rite - 1].value > peak.value*1e-2)
				break;

		Quantity iPC = Math2.quadargmax(left, rite, this.arealDensity); // index of max compression
		Quantity peakCompress = Math2.interp(timeAxis, iPC); // time of max compression
		Quantity iTPeak = Math2.quadargmax(left, rite, this.ionTemperature);
		Quantity[] moments = new Quantity[5];
		for (int k = 0; k < moments.length; k ++)
			moments[k] = Math2.moment(k, timeBins, this.neutronYield);
		Quantity burnWidth = Math2.fwhm(timeBins, this.neutronYield);

		Quantity[] res = {
			  new Quantity((endTime - startTime)/1000., dofs),
			  moments[0].times(timeStep), bangTime,
			  burnWidth, moments[3], moments[4],
			  peakCompress.minus(bangTime),
			  Math2.average(this.ionTemperature, this.neutronYield, left, rite),
			  Math2.quadInterp(this.ionTemperature, iTPeak),
			  Math2.quadInterp(this.ionTemperature, iPC),
			  Math2.quadInterp(this.ionTemperature, iBT),
			  Math2.derivative(timeAxis, this.ionTemperature, iPC, .10, 1),
			  Math2.derivative(timeAxis, this.ionTemperature, bangTime, .10, 1),
			  Math2.derivative(timeAxis, this.ionTemperature, bangTime, .10, 2),
			  Math2.average(this.arealDensity, this.neutronYield, left, rite),
			  Math2.quadInterp(this.arealDensity, iPC),
			  Math2.quadInterp(this.arealDensity, iBT),
			  Math2.derivative(timeAxis, this.arealDensity, bangTime, .10, 1),
			  Math2.derivative(timeAxis, V, bangTime, .12, 2).over(Math2.quadInterp(V, iBT)),
		}; // collect the figures of merit

		logger.info(String.format("Total yield (μ0):  %s", res[1].times(1e15).toString(covarianceMatrix)));
		logger.info(String.format("Bang time:         %s ns", res[2].toString(covarianceMatrix)));
		logger.info(String.format("Burn width (μ2):   %s ps", res[3].over(1e-3).toString(covarianceMatrix)));
		logger.info(String.format("Burn skewness (μ3):%s", res[4].toString(covarianceMatrix)));
		logger.info(String.format("Burn kurtosis (μ4):%s", res[5].toString(covarianceMatrix)));
		logger.info(String.format("Peak compression:  %s ps + BT", res[6].over(1e-3).toString(covarianceMatrix)));
		logger.info(String.format("Burn-averaged Ti:  %s keV", res[7].toString(covarianceMatrix)));
		logger.info(String.format("Ti at peak:        %s keV", res[8].toString(covarianceMatrix)));
		logger.info(String.format("Ti at compression: %s keV", res[9].toString(covarianceMatrix)));
		logger.info(String.format("Ti at BT:          %s keV", res[10].toString(covarianceMatrix)));
		logger.info(String.format("dTi/dt at compr.:  %s g/cm^2/(100 ps)", res[11].over(1e1).toString(covarianceMatrix)));
		logger.info(String.format("dTi/dt at BT:      %s keV/(100 ps)", res[12].over(1e1).toString(covarianceMatrix)));
		logger.info(String.format("d^2Ti/dt^2 at BT:  %s keV/ns^2", res[13].toString(covarianceMatrix)));
		logger.info(String.format("Burn-averaged \u03C1R:  %s g/cm^2", res[14].toString(covarianceMatrix)));
		logger.info(String.format("\u03C1R at compression: %s g/cm^2", res[15].toString(covarianceMatrix)));
		logger.info(String.format("\u03C1R at BT:          %s g/cm^2", res[16].toString(covarianceMatrix)));
		logger.info(String.format("d\u03C1R/dt at BT:      %s g/cm^2/(100 ps)", res[17].over(1e1).toString(covarianceMatrix)));
		logger.info(String.format("d^2V/dt^2/V at BT: %s 1/ns^2", res[18].toString(covarianceMatrix)));

		double[] output = new double[2*res.length]; // and error bars and serve
		for (int i = 0; i < res.length; i++) {
			output[2*i + 0] = res[i].value;
			output[2*i + 1] = Math.sqrt(res[i].variance(covarianceMatrix));
		}
		assert output.length == HEADERS_WITH_ERRORS.length - 1;
		return output;
	}

	/**
	 * take the time-corrected spectrum and use it to compute and store time-
	 * resolved values for ion temperature, areal density, and yield.  also store
	 * the covariance matrix for them.
	 * @param signal the time-corrected spectrum we want to understand
	 * @param errorMode should error bars be computed (it's kind of intensive)?
	 */
	private void analyze(double[][] signal, ErrorMode errorMode) {
		for (double[] doubles: signal)
			for (double value: doubles)
				if (Double.isNaN(value))
					throw new IllegalArgumentException("well fuck you too");

		final int N = 3;
		final int M = timeAxis.length;

		// start with the simplest possible reconstruction
		double[] naiveNeutronYield = new double[signal[0].length];
		double efficiency = this.averageEfficiency(4, .1, 0, Double.POSITIVE_INFINITY);
		for (int i = 0; i < signal.length; i ++) {
			double background = this.detector.background(energyAxis[i], energyBins, timeBins, true);
			for (int j = 0; j < signal[i].length; j ++) {
				naiveNeutronYield[j] += Math.max(
						0, (signal[i][j] - background)/detector.gain/efficiency)/(timeStep*1e15); // the spectrum scaled by the efficiency of the system
			}
		}

		double[] yieldScale = new double[M];
		for (int j = 0; j < M; j++)
			yieldScale[j] = Math2.max(naiveNeutronYield)/6.;
		double[] temperatureScale = Math2.full(10, M);
		double[] densityScale = Math2.full(1, M);
		double[] lowerBound = Math2.full(0, M);
		double[] upperBound = Math2.full(Double.POSITIVE_INFINITY, M);

		// put together the initial gess for the proper fitting
		double[] ionTemperature = Math2.full(4, M); // keV
		double[] electronTemperature = Math2.full(ELECTRON_TEMPERATURE, M); // keV
		double[] bulkFlowVelocity = Math2.full(BULK_FLOW_VELOCITY, M); // μm/ns
		double[] arealDensity = Math2.full(100e-3, M); // g/cm^2

		// then do the second simplest reconstruction for the burn history
		double[] finalIonTemperature = ionTemperature;
		double[] finalArealDensity = arealDensity;
		double[] neutronYield = optimize(
			  (double[] x) -> logPosterior(
					signal,
					x,
					finalIonTemperature,
					electronTemperature,
					bulkFlowVelocity,
					finalArealDensity,
					null, null, true, 100),
			  naiveNeutronYield, yieldScale, lowerBound, upperBound); // 10^15/ns

//				try {
//					double[][] trueTrajectories = CSV.read(new File("input/trajectories og with falling temp.csv"), ',', 1);
//					double[] time = new double[trueTrajectories.length];
//					double[] ρR = new double[trueTrajectories.length];
//					double[] Yn = new double[trueTrajectories.length];
//					double[] Ti = new double[trueTrajectories.length];
//					for (int i = 0; i < trueTrajectories.length; i++) {
//						time[i] = trueTrajectories[i][0];
//						Yn[i] = 10*trueTrajectories[i][1]*.1*1e6/1e-6/(14e6*1.6e-19)/1e15*1e-9; // convert from 0.1MJ/μs to 1e15n/ns
//						Ti[i] = trueTrajectories[i][4];
//						ρR[i] = trueTrajectories[i][3];
//					}
//					for (int j = 0; j < timeAxis.length; j ++) {
//						if (timeAxis[j] < time[0]) {
//							neutronYield[j] = 0;
//							ionTemperature[j] = Ti[0];
//							arealDensity[j] = ρR[0];
//						} else if (timeAxis[j] > time[time.length-1]) {
//							neutronYield[j] = 0;
//							ionTemperature[j] = Ti[time.length-1];
//							arealDensity[j] = ρR[time.length-1];
//						} else {
//							neutronYield[j] = Math2.interp(timeAxis[j], time, Yn);
//							ionTemperature[j] = Math2.interp(timeAxis[j], time, Ti);
//							arealDensity[j] = Math2.interp(timeAxis[j], time, ρR);
//						}
//					}
//					System.out.println("took the correct anser as the initial gess");
//				}
//				catch (IOException e) {
//					System.err.println("the test thing didn't work.");
//				}

		// determine the bounds of the region that's worth optimizing
		double[] statistics = new double[M];
		for (int j = 0; j < M; j ++)
			statistics[j] = neutronYield[j]*1e15*efficiency*timeStep;
		final int peak = Math.max(1, Math.min(statistics.length - 2,
		                                      Math2.argmax(statistics)));
		int left = -1;
		for (int j = peak - 2; j >= -1; j --) { // scan backward for the start of the signal
			if (left == -1 && (j == -1 || statistics[j] < 10)) // if the statistics are very lo
				left = j + 1; // this is probably the start
			else if (j >= 0 && statistics[j] > 100) // if they earlier go somewhat hi, tho,
				left = -1; // then we were mistaken
		}
		int rite = -1;
		for (int j = peak + 2; j <= statistics.length; j ++) { // scan forward for the end of the signal
			if (rite == -1 && (j == statistics.length || statistics[j] < 10)) // if the statistics are very lo
				rite = j; // this is probably the end
			else if (j < statistics.length && statistics[j] > 100) // if they later go somewhat hi, tho,
				rite = -1; // then we were mistaken
		}
		boolean[] active = new boolean[timeAxis.length];
		for (int j = 0; j < timeAxis.length; j ++) {
			active[j] = j >= left && j < rite;
			if (!active[j]) neutronYield[j] = 0;
		}
		System.out.println("active is");
		System.out.println(Arrays.toString(active));
		System.out.println("because of");
		System.out.println(Arrays.toString(statistics));

		final double finalSmoothing = 1;//smoothing;

		boolean[][] signalSpread = Math2.nonzero(
				this.detector.response(
						energyBins, timeBins, this.ionOptics.response(
								energyBins, timeBins, SpectrumGenerator.generateSpectrum(
										neutronYield, ionTemperature, electronTemperature, bulkFlowVelocity, arealDensity,
										energyBins, timeBins),
								false, false),
						false, true, true));

		double lastValue;
		double thisValue = logPosterior(
			  signal, neutronYield,
			  ionTemperature, electronTemperature,
			  bulkFlowVelocity, arealDensity,
			  active, signalSpread, false, finalSmoothing);

		do {
			final double[] neutronYieldInitialGess = neutronYield;
			final double[] ionTemperatureInitialGess = ionTemperature;
			final double[] arealDensityInitialGess = arealDensity;

			logger.log(Level.FINE, "Fitting ion temperature...");

			final double[] ionTemperatureNewGess = optimize(
				  (double[] x) -> logPosterior(
						signal,
						neutronYieldInitialGess,
						x,
						electronTemperature,
						bulkFlowVelocity,
						arealDensityInitialGess,
						active, signalSpread, false, finalSmoothing),
				  ionTemperatureInitialGess,
				  temperatureScale, lowerBound, upperBound, active);

			logger.log(Level.FINE, "Fitting ρR trajectory...");

			final double[] arealDensityNewGess = optimize(
				  (double[] x) -> logPosterior(
						signal,
						neutronYieldInitialGess,
						ionTemperatureNewGess,
						electronTemperature,
						bulkFlowVelocity,
						x,
						active, signalSpread, false, finalSmoothing),
				  arealDensityInitialGess,
				  densityScale, lowerBound, upperBound, active);

			logger.log(Level.FINE, "Fitting burn history...");

			final double[] neutronYieldNewGess = optimize(
				  (double[] x) -> logPosterior(
						signal,
						x,
						ionTemperatureNewGess,
						electronTemperature,
						bulkFlowVelocity,
						arealDensityNewGess,
						active, signalSpread, false, finalSmoothing),
				  neutronYieldInitialGess,
				  yieldScale, lowerBound, upperBound, active);

			ionTemperature = ionTemperatureNewGess;
			arealDensity = arealDensityNewGess;
			neutronYield = neutronYieldNewGess;

			lastValue = thisValue;
			thisValue = logPosterior(
				  signal, neutronYield,
				  ionTemperature, electronTemperature,
				  bulkFlowVelocity, arealDensity,
				  active, signalSpread,false, finalSmoothing);
		} while (lastValue - thisValue > this.precision);

		this.neutronYield = new Quantity[M];
		this.ionTemperature = new Quantity[M];
		this.arealDensity = new Quantity[M];
		for (int j = 0; j < timeAxis.length; j ++) {
			this.neutronYield[j] = new Quantity(
					neutronYield[j], j, N*M);
			this.ionTemperature[j] = new Quantity(
					ionTemperature[j], M+j, N*M);
			this.arealDensity[j] = new Quantity(
					arealDensity[j], 2*M+j, N*M);
		}

		if (errorMode == ErrorMode.HESSIAN) { // calculate the error bars using second-derivatives
			logger.log(Level.INFO, "calculating error bars.");
			final int j0 = left, m = rite - left;
			double[] x0 = new double[N*m], dx = new double[N*m];
			for (int j = left; j < rite; j ++) { // set up for a finite-difference hessian calculation
				x0[      j - j0] = neutronYield[j];
				dx[      j - j0] = Math2.max(neutronYield)*1e-4;
				x0[  m + j - j0] = ionTemperature[j];
				dx[  m + j - j0] = 1e-2;
				x0[2*m + j - j0] = arealDensity[j];
				dx[2*m + j - j0] = 1e-3;
			}
			final double[] testNeutronYield = Arrays.copyOf(neutronYield, M);
			final double[] testIonTemperature = Arrays.copyOf(ionTemperature, M);
			final double[] testArealDensity = Arrays.copyOf(arealDensity, M);
			double[][] hessian = Math2.hessian((double[] x) -> {
				assert x.length == N*m;
				System.arraycopy(x, 0  , testNeutronYield, j0, m);
				System.arraycopy(x, 1*m, testIonTemperature, j0, m);
				System.arraycopy(x, 2*m, testArealDensity, j0, m);
				return this.logPosterior(signal, testNeutronYield,
										 testIonTemperature, electronTemperature,
										 bulkFlowVelocity, testArealDensity,
										 active, signalSpread, false, finalSmoothing);
			}, x0, dx); // the do that
			Math2.coercePositiveSemidefinite(hessian);
			double[][] activeCovarianceMatrix = Math2.pseudoinv(hessian); // invert it to get the covariance
			for (int i = 0; i < hessian.length; i ++) {
				if (activeCovarianceMatrix[i][i] < 1/hessian[i][i]) // these are all approximations, and sometimes they violate the properties of positive semidefiniteness
					activeCovarianceMatrix[i][i] = 1/hessian[i][i]; // do what you must to make it work
			}
			covarianceMatrix = new double[N*M][N*M]; // then expand it to include the inactive dimensions
			for (int k1 = 0; k1 < N; k1 ++)
				for (int j1 = left; j1 < rite; j1++)
					for (int k2 = 0; k2 < N; k2++)
						for (int j2 = left; j2 < rite; j2++)
							covarianceMatrix[k1*M + j1][k2*M + j2] = activeCovarianceMatrix[k1*m + j1 - j0][k2*m + j2 - j0];
		}
		else if (errorMode == ErrorMode.STATISTICS) {
			covarianceMatrix = new double[N*M][N*M];
			for (int j = 0; j < rite - left; j ++) {
				double primaryStatistics = 1 +
						neutronYield[j]*1e15*
								ionOptics.efficiency(14)*
								timeStep; // total signal deuteron yield from this neutron time bin
				double dsStatistics = 1 + primaryStatistics*arealDensity[j]/21.;
				covarianceMatrix[    j][    j] = Math.pow(neutronYield[j], 2)/primaryStatistics;
				covarianceMatrix[  M+j][  M+j] = Math.pow(ionTemperature[j], 2)*2/(primaryStatistics - 1);
				covarianceMatrix[2*M+j][2*M+j] = Math.pow(arealDensity[j], 2)/dsStatistics;
			}
		}
		else {
			covarianceMatrix = new double[N*M][N*M];
		}
	}

	/**
	 * a utility function for the fitting: get the error in this spectrum fit
	 * @param active whether to apply the prior to each timestep
	 * @param sumInEnergy whether to combine energies when comparing
	 * @param signalRegion only pay attention to pixels markd by signalRegion
	 * @return the log of the inverse posterior value for this spectrum
	 */
	private double logPosterior(double[][] spectrum,
								double[] neutronYield,
								double[] ionTemperature,
								double[] electronTemperature,
								double[] bulkFlowVelocity,
								double[] arealDensity,
								boolean[] active,
								boolean[][] signalRegion,
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
		double[][] deuterons = this.ionOptics.response(
			  energyBins, timeBins, neutrons,
			  false, false);
		double[][] signal = this.detector.response(
			  energyBins, timeBins, deuterons,
			  false, true, true);

		double totalError = 0; // negative log likelihood
		for (int j = 0; j < timeAxis.length; j ++) {
			int numEnergies = (sumInEnergy) ? 1 : spectrum.length;
			double[] experValues = new double[numEnergies];
			double[] theorValues = new double[numEnergies];
			double[] variances = new double[numEnergies];
			double[] backgrounds = new double[numEnergies];
			for (int i = 0; i < spectrum.length; i ++) {
				if (signalRegion == null || signalRegion[i][j]) {
					int index = (sumInEnergy) ? 0 : i;
					experValues[index] += spectrum[i][j];
					theorValues[index] += signal[i][j];
					variances[index] += detector.noise(energyAxis[i], energyBins, timeBins, true);
					backgrounds[index] = detector.background(energyAxis[i], energyBins, timeBins, true);
				}
			}
			for (int i = 0; i < numEnergies; i ++) {
				double theorNumber = (theorValues[i] - backgrounds[i])/detector.gain;
				double experNumber = (experValues[i] - backgrounds[i])/detector.gain;
				if (variances[i] > 0) { // if this detector has significant noise
					double variance = variances[i]
							+ theorNumber*detector.gain // include the poisson noise of the signal electrons
							+ theorNumber*detector.gain*detector.gain; // include pre-amplification poisson noise
					totalError += Math.pow(experValues[i] - theorValues[i], 2)/
							(2*variance); // and use a Gaussian approximation
				}
				else { // if the detector noise is zero
					assert !Double.isNaN(theorNumber);
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
		for (int j = -2; j < timeAxis.length + 2; j ++) {
			double[] y = new double[3];
			for (int k = 0; k < 3; k ++) {
				if (j+k >= 0 && j+k < timeAxis.length && (active == null || active[j+k]))
					y[k] = neutronYield[j+k];
				else
					y[k] = 0;
			}
			double slope0 = (y[0] != 0 || y[1] != 0) ? (y[1] - y[0])/(y[0] + y[1]) : 0;
			double slope1 = (y[1] != 0 || y[2] != 0) ? (y[2] - y[1])/(y[1] + y[2]) : 0;
			totalPenalty += smoothing*1e-0/timeStep*
					(Math.exp(slope1 - slope0) - Math.exp((slope1 - slope0)/2)*2 + 1); // encourage a smooth burn history with no local mins
		}
		double totalYield = Math2.sum(neutronYield)*timeStep;
		for (int j = 0; j < timeAxis.length; j ++)
			totalPenalty -= 0e-1*Math.pow(neutronYield[j]/totalYield, 2)*timeStep; // encourage a peaked Yn (short burn width)
		for (int j = 0; j < timeAxis.length; j ++)
			totalPenalty += arealDensity[j]/5.0 - Math.log(arealDensity[j])/50.0; // gamma prior on areal density
		double expected_temperature = 4;//5.5e-4*Math.pow(totalYield*1e15, .25);
		for (int j = 0; j < timeAxis.length; j ++)
			totalPenalty += Math.pow((Math.log(ionTemperature[j]/expected_temperature))/2.5, 2); // log-normal prior on Ti
//		double expected_temperature_slope = 6.0*(Math.log10(totalYield) + 15) - 100.5;
		for (double[] x: new double[][] {ionTemperature, arealDensity}) {
			for (int j = 0; j < timeAxis.length - 2; j ++) {
				if (active == null || (active[j] && active[j+1] && active[j+2])) {
					double Ψpp = (x[j] - 2*x[j+1] + x[j+2])/
							Math.pow(timeStep, 3);
					double Ψ = (x[j] + x[j+1] + x[j+2])/3;
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

	public double[] getDeuteronEnergyBins() {
		return this.deuteronEnergyBins;
	}

	public double[][] getDeuteronSpectrum() {
		return this.deuteronSpectrum;
	}

	public double[][] getSignalDistribution() {
		return this.signalDistribution;
	}

	public double[][] getFitNeutronSpectrum() {
		return this.fitNeutronSpectrum;
	}

	public double[][] getFitDeuteronSpectrum() {
		return this.fitDeuteronSpectrum;
	}

	public double[][] getFitSignalDistribution() {
		return this.fitSignalDistribution;
	}

	public double[][] getIdealSignalDistribution() {
		return this.idealSignalDistribution;
	}

	public double[][] getBackgroundSpectrum() {
		double[][] background = new double[energyAxis.length][timeAxis.length];
		for (int i = 0; i < energyAxis.length; i ++)
			background[i] = Math2.full(
				  detector.background(energyAxis[i], energyBins, timeBins, true),
				  timeAxis.length);
		return background;
	}

	public double[][] efficiencyCorrect(double[][] distribution, boolean gaps) {
		double[][] spectrum = new double[distribution.length][distribution[0].length];
		for (int i = 0; i < spectrum.length; i ++) {
			for (int j = 0; j < spectrum[i].length; j ++) {
				double z = distribution[i][j];
				double z0 = detector.background(energyAxis[i], energyBins, timeBins, gaps);
				double η = detector.efficiency(energyAxis[i], gaps)*detector.gain;
				spectrum[i][j] = (z - z0)/η;
			}
		}
		return spectrum;
	}

	/**
	 * get the time bin centers in ns
	 */
	public double[] getTimeAxis() {
		return this.timeAxis;
	}

	/**
	 * get the energy bin centers in MeV
	 */
	public double[] getEnergyAxis() {
		return this.energyAxis;
	}

	/**
	 * get the energy bin centers in MeV
	 */
	public double[] getDeuteronEnergyAxis() {
		return this.deuteronEnergyAxis;
	}

	/**
	 * get the neutron yield mean values in [10^15/ns]
	 */
	public double[] getNeutronYield() {
		return Math2.modes(this.neutronYield);
	}

	/**
	 * get the neutron yield error bars in [10^15/ns]
	 */
	public double[] getNeutronYieldError() {
		return Math2.stds(this.neutronYield, this.covarianceMatrix);
	}

	/**
	 * get the ion temperature mean values in [keV]
	 */
	public double[] getIonTemperature() {
		return Math2.modes(this.ionTemperature);
	}

	/**
	 * get the ion temperature error bars in [keV]
	 */
	public double[] getIonTemperatureError() {
		return Math2.stds(this.ionTemperature, this.covarianceMatrix);
	}

	/**
	 * get the ρR mean values in [g/cm^2]
	 */
	public double[] getArealDensity() {
		return Math2.modes(this.arealDensity);
	}

	/**
	 * get the ρR error bars in [g/cm^2]
	 */
	public double[] getArealDensityError() {
		return Math2.stds(this.arealDensity, this.covarianceMatrix);
	}


	private static String[] appendErrorsToHeader() {
		String[] out = new String[2*HEADERS.length - 1];
		for (int i = 0; i < HEADERS.length; i ++) {
			if (i < 2)
				out[i] = HEADERS[i];
			else {
				out[2*(i-1)+1] = HEADERS[i];
				out[2*(i-1)+2] = HEADERS[i] + " error";
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
