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

import util.COSYMapping;
import util.CSV;
import util.Math2;
import util.Math2.DiscreteFunction;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;

import static physics.Analysis.MC_RANDOM;
import static physics.Analysis.NOISE_RANDOM;

/**
 * A class to handle all of the deuteron physics, from their birth in the foil to
 * their collision with the detector plane.
 */
public class IonOptics {

	public static class IonOpticConfiguration {
		public static IonOpticConfiguration HIGH_EFFICIENCY =
			  new IonOpticConfiguration(90e-6, 400e-6, 5e-3);
		public static IonOpticConfiguration MID_EFFICIENCY =
			  new IonOpticConfiguration(40e-6, 200e-6, 3e-3);
		public static IonOpticConfiguration LOW_EFFICIENCY =
			  new IonOpticConfiguration(25e-6, 100e-6, 2e-3);

		public final double foilThickness;
		public final double foilRadius;
		public final double apertureWidth;

		public IonOpticConfiguration(
			  double foilThickness, double foilRadius, double apertureWidth) {
			this.foilThickness = foilThickness;
			this.foilRadius = foilRadius;
			this.apertureWidth = apertureWidth;
		}
	}

	private static final double keV = 1e3*Math.abs(Particle.E.charge);
	private static final double MeV = 1e6*Math.abs(Particle.E.charge);
	private static final double μm = 1e-6;
	private static final double ns = 1e-9;

	private static final int x = 0, y = 1, z = 2;

	private static final int TIME_CORRECTION_RESOLUTION = 50;
	private static final int STOPPING_DISTANCE_RESOLUTION = 64;
	private static final String STOPPING_POWER_FILENAME = "input/stopping_power_%ss_CD.csv";

	private static final int TRANSFER_MATRIX_TRIES = 10000; // the number of points to sample in each column of the transfer matrix

	private final double foilDistance; // z coordinate of midplane of foil [m]
	private final double foilWidth; // extent of foil in dispersion direccion [m]
	private final double foilHeight; // extent of foil in nondispersion direccion [m]
	private final double foilThickness; // thickness of foil [m]
	private final double apertureDistance; // distance from TCC to aperture [m]
	private final double apertureWidth; // horizontal dimension of aperture [m]
	private final double apertureHeight; // vertical dimension of aperture [m]
	private final double minEd, maxEd; // bounds on deuterons that will be detected by the CSI [J]
	private final double focalPlaneAngle; // angle of focal plane [rad]
	private final COSYMapping cosyMapping; // table of polynomials generated by COSY
	private final double calibrationPrecision; // the error in the transfer function
	private final boolean reuseMatrix; // if this is true, it ignores all this and loads whatever it had last time

	public final double energyFactor; // conversion factor between neutron and ion energies []
	private final double probHitsFoil; // probability that the neutron goes through the foil

	private double[] energyBins; // the energy bins that were last used to compute the transfer matrix
	private double[] timeBins; // the time bins that were last used to compute the transfer matrix
	private double[][] rongTransferMatrix; // the theoretical transfer matrix of the full ion-optick system
	private double[][] trueTransferMatrix; // a randomly varied transfer matrix, to mimick the uncertainty of reality

	private final DiscreteFunction distanceVsEnergy; // stopping distance info
	private final DiscreteFunction energyVsDistance; // inverse stopping distance info
	private final DiscreteFunction energyVsPosition; // map between location on detector and energy going into lens

	/**
	 * generate an IonOptic object from a configuration preset
	 * @param config the foil/aperture configuration
	 * @param tiltAngle the angle of the detector (degrees)
	 * @param precision how accurately we expect to know the resolution (ratio)
	 * @param reuseMatrix whether to load the nominal response matrix from disk
	 */
	public IonOptics(IonOpticConfiguration config,
					 String cosyFile, double tiltAngle,
					 double precision, boolean reuseMatrix) throws IOException {
		this(3.0e-3,
			 2*config.foilRadius,
			 2*config.foilRadius,
			 config.foilThickness,
			 6.0e+0,
			 config.apertureWidth,
			 20.0e-3,
			 12, 16,
			 CSV.readCosyCoefficients(new File(String.format(
				   "input/%s.txt", cosyFile)),
				  3, Particle.D, 12.45),
			 tiltAngle,
			 precision,
			 reuseMatrix);
	}

	/**
	 * put together the ion optic simulacion
	 * @param foilDistance the distance from TCC to the foil (m)
	 * @param foilWidth the total width of the foil (m)
	 * @param foilHeight the total hite of the foil (m)
	 * @param foilThickness the thickness of the foil (m)
	 * @param apertureDistance the distance from TCC to the aperture (m)
	 * @param apertureWidth the width of the aperture (m)
	 * @param apertureHeight the hite of the aperture (m)
	 * @param minEn the loest neutron energy to hit the detector (MeV)
	 * @param maxEn the hiest neutron energy to hit the detector (MeV)
	 * @param cosyMapping the COSY coefficient matrix
	 * @param focalTilt angle of the focal plane (0 means untilted) (deg)
	 * @throws IOException if it can't find the stopping power file
	 * @throws NumberFormatException if the stopping power file accepts bribes
	 */
	public IonOptics(
			double foilDistance, double foilWidth, double foilHeight, double foilThickness,
	        double apertureDistance, double apertureWidth, double apertureHeight,
	        double minEn, double maxEn, COSYMapping cosyMapping, double focalTilt,
			double calibrationPrecision, boolean reuseMatrix) throws IOException {

		Particle ion = cosyMapping.ion;
		double A = ion.mass/Particle.N.mass;
		this.energyFactor = 4*A/Math.pow(A + 1, 2);
		this.foilDistance = foilDistance;
		this.foilWidth = foilWidth;
		this.foilHeight = foilHeight;
		this.foilThickness = foilThickness;
		this.apertureDistance = apertureDistance;
		this.apertureWidth = apertureWidth;
		this.apertureHeight = apertureHeight;
		this.minEd = energyFactor*minEn*MeV;
		this.maxEd = energyFactor*maxEn*MeV;
		this.focalPlaneAngle = Math.toRadians(focalTilt);
		this.cosyMapping = cosyMapping;
		this.probHitsFoil = foilWidth*foilHeight/(4*Math.PI*foilDistance*foilDistance);
		this.calibrationPrecision = calibrationPrecision;
		this.reuseMatrix = reuseMatrix;

		double[][] stoppingData = CSV.read(
				new File(String.format(STOPPING_POWER_FILENAME, ion.name)),
				',');
		for (int i = 0; i < stoppingData.length; i ++) {
			stoppingData[i][1] = 1/(stoppingData[i][1]*keV/μm); // converting from [keV/μm]
			stoppingData[i][0] = stoppingData[i][0]*keV; // and from [keV]
		}
		DiscreteFunction distanceVsEnergyRaw = new DiscreteFunction(stoppingData, true).antiderivative(); // integrate the stopping power to get stopping distance
		this.distanceVsEnergy = distanceVsEnergyRaw.indexed(STOPPING_DISTANCE_RESOLUTION); // m(J)
		this.energyVsDistance = distanceVsEnergyRaw.inv().indexed(STOPPING_DISTANCE_RESOLUTION); // J(m)

		double[] calibEnergies = new double[2*TIME_CORRECTION_RESOLUTION];
		double[] detectorPosition = new double[calibEnergies.length];
		for (int i = 0; i < calibEnergies.length; i ++) {
			calibEnergies[i] = (minEd + (this.maxEd - minEd)*i/(calibEnergies.length-1));
			double[] v = {0, 0, Math.sqrt(2*calibEnergies[i]/ion.mass)};
			double[] r = computeFocusedPosition(new double[] {0,0,0}, v, 0);
			detectorPosition[i] = r[x]/Math.cos(focalPlaneAngle);
		}
		this.energyVsPosition = new DiscreteFunction(calibEnergies, detectorPosition).inv()
				.indexed(TIME_CORRECTION_RESOLUTION); // (m -> J deuteron)
	}

	/**
	 * compute the probability that a given neutron released at this energy will spawn an ion
	 * and knock that ion through the aperture.
	 * @param energy energy of the released particles [MeV]
	 * @return the fraction of particles that are worth simulating
	 */
	public double efficiency(double energy) {
		double n = 0.08e2; // I'm not sure what units this has or whence it came
		double dσdΩ = 4.3228/Math.sqrt(energy) - 0.6523; // same with these ones
		double dΩ = apertureWidth*apertureHeight / Math.pow(apertureDistance - foilDistance, 2);
		return probHitsFoil * n*dσdΩ*dΩ*foilThickness; // assume the foil is thin so we don't have to worry about multiple collisions
	}

	/**
	 * compute the dispersion at the detector plane
	 * @return the time skew [keV/mm]
	 */
	public double computeDispersion() {
		double[] central = this.backCalculate(0, 0, 0, 0);
		double[] perturbd = this.backCalculate(
			  1e-6*Math.cos(focalPlaneAngle), 1e-6*Math.sin(focalPlaneAngle), 0, 0);
		return (perturbd[0] - central[0])/1e-6;
	}

	/**
	 * compute the time skew at the detector plane
	 * @return the time skew [ps/keV]
	 */
	public double computeTimeSkew() {
		double[] central = this.backCalculate(0, 0, 0, 0);
		double[] perturbd = this.backCalculate(
			  1e-6*Math.cos(focalPlaneAngle), 1e-6*Math.sin(focalPlaneAngle), 0, 0);
		return (perturbd[1] - central[1])/1e-12/((perturbd[0] - central[0])*1e3);
	}

	/**
	 * compute the time and energy resolution (FWHM) for particles at a given energy.
	 * @param referenceEnergy the neutron energy of the central ray [MeV]
	 * @return {energy resolution [keV], time resolution [ps]}
	 */
	public double[] computeResolution(double referenceEnergy) {
		int nEnergies = 30;
		int nTimes = 20;
		double[] energyRange = {-1.1, 0.1}; // [MeV]
		double[] timeRange = {-.160, .160}; // [ns]

		double[] energyBins = new double[nEnergies+1]; // [MeV]
		for (int i = 0; i <= nEnergies; i ++)
			energyBins[i] = (energyRange[1] - energyRange[0])*i/nEnergies +
					energyRange[0] + referenceEnergy;
		double[] timeBins = new double[nTimes+1]; // [ns]
		for (int j = 0; j <= nTimes; j ++)
			timeBins[j] = (timeRange[1] - timeRange[0])*j/nTimes + timeRange[0];

		makeSureWeHaveTransferMatrix(energyBins, timeBins);

		int iRef = Math2.bin(referenceEnergy, energyBins);
		int jRef = nTimes/2;

		double[] energyDist = new double[nEnergies];
		double[] timeDist = new double[nTimes];
		for (int i = 0; i < nEnergies; i ++) {
			for (int j = 0; j < nTimes; j ++) {
				energyDist[i] += rongTransferMatrix[(timeBins.length-1)*i + j][(timeBins.length-1)*iRef + jRef];
				timeDist[j] += rongTransferMatrix[(timeBins.length-1)*i + j][(timeBins.length-1)*iRef + jRef];
			}
		}

		//		double[] energyAxis = new double[nEnergies];
		//		for (int i = 0; i < nEnergies; i ++)
		//			energyAxis[i] = (energyBins[i] + energyBins[i+1])/2.;
		//		double[] timeAxis = new double[nTimes];
		//		for (int j = 0; j < nTimes; j ++)
		//			timeAxis[j] = (timeBins[j] + timeBins[j+1])/2.;

		//		System.out.println(Arrays.toString(energyAxis));
//		System.out.println(Arrays.toString(energyDist));
//		System.out.println(Arrays.toString(timeAxis));
//		System.out.println(Arrays.toString(timeDist));

		double energyResolution = Math2.std(energyBins, energyDist);
		double timeResolution = Math2.std(timeBins, timeDist);
		return new double[] { energyResolution/1e-3, timeResolution/1e-3};
	}

	/**
	 * compute the total response to an implosion with the given neutron spectrum using the
	 * precomputed transfer matrix. account for electrostatic time correction, but not for any
	 * analysis.
	 * @param energyBins the edges of the energy bins [MeV]
	 * @param timeBins the edges of the time bins [ns]
	 * @param inSpectrum the neutron counts in each bin
	 * @param stochastic whether to add noise to mimic real data
	 * @param actual whether to use the true response function or what we believe to be the response function
	 * @return the counts in the measured spectrum bins
	 */
	public double[][] response(double[] energyBins, double[] timeBins, double[][] inSpectrum,
	                           boolean stochastic, boolean actual) {
		makeSureWeHaveTransferMatrix(energyBins, timeBins); // the full nmxnm believed transfer matrix

		double[] u = new double[(energyBins.length-1)*(timeBins.length-1)];
		for (int i = 0; i < energyBins.length-1; i ++)
			for (int j = 0; j < timeBins.length-1; j ++)
				u[(timeBins.length-1)*i + j] = inSpectrum[i][j]; // now flatten the spectrum

		double[] v;
		if (actual)
			v = Math2.matmul(trueTransferMatrix, u);
		else
			v = Math2.matmul(rongTransferMatrix, u); // then do the multiplication

		double[][] outSpectrum = new double[energyBins.length-1][timeBins.length-1];
		for (int i = 0; i < energyBins.length-1; i ++)
			for (int j = 0; j < timeBins.length-1; j ++)
				outSpectrum[i][j] = v[(timeBins.length-1)*i + j]; // finally, unflatten. inflate. rounden.

		if (stochastic) { // to simulate stochasticity
			for (int i = 0; i < energyBins.length-1; i ++)
				for (int j = 0; j < timeBins.length-1; j ++)
					outSpectrum[i][j] = Math2.poisson(outSpectrum[i][j], NOISE_RANDOM); // just jitter every cell
		}

		return outSpectrum;
	}

	/**
	 * simulate a single random neutron emitted from TCC at the given energy and determine the
	 * position and time at which its child ion crosses the focal plane, time corrected, and
	 * the back-calculated energy.
	 * @param energy initial energy of released neutron [MeV].
	 * @param time initial time of released neutron [s].
	 * @return { x, y, z, t } [m, m, m, s].
	 */
	public double[] simulate(double energy, double time, boolean foilBlur) {
		double[] rCollision = chooseCollisionPosition();

		double[] rAperture = chooseAperturePosition();

		double[] vFinal = computeFinalVelocity(energy, rCollision, rAperture, foilBlur);

		return computeFocusedPosition(rCollision, vFinal, time);
	}

	/**
	 * simulate a single neutron emitted from TCC at the given energy and determine the
	 * position and time at which its child ion crosses the focal plane, time corrected, and
	 * the back-calculated energy.
	 * @param energy initial energy of released neutron (MeV).
	 * @return { x, y, z, t } [m, m, m, s].
	 */
	public double[] map(double energy) {
		energy = energyFactor*energy*MeV;
		double[] r = {0, 0, 0};
		double[] vFinal = {0, 0, Math.sqrt(2*energy/cosyMapping.ion.mass)};
		return computeFocusedPosition(r, vFinal, 0);
	}

	/**
	 * evaluate the matrix governing the response of the MRSt to neutrons at
	 * particular energies <i>if</i> we haven't yet for these bins.  if we have
	 * already evaluated it, do absolutely noting.
	 * @param energyBins the neutron energy bins by which we linearize [MeV]
	 * @param timeBins the neutron time bins by which we linearize [ns]
	 */
	private void makeSureWeHaveTransferMatrix(double[] energyBins, double[] timeBins) {
		if (this.rongTransferMatrix == null ||
			  !Arrays.equals(energyBins, this.energyBins) ||
			  !Arrays.equals(timeBins, this.timeBins)) { // (only evaluate it if it hasn't been evaluated yet)
			this.energyBins = energyBins;
			this.timeBins = timeBins;

			this.rongTransferMatrix = null;

			if (this.reuseMatrix) { // if the user sed to reuse the last transfer matrix
				try {
					this.rongTransferMatrix = CSV.read(
						  new File("output/transfer matrix.csv"), ','); // load it from disc
				} catch (IOException e) {
					e.printStackTrace();
				}
				int n = (energyBins.length - 1)*(timeBins.length - 1);
				if (this.rongTransferMatrix.length != n || this.rongTransferMatrix[0].length != n) { // check to make sure it is correctly shaped
					System.out.println("I couldn't reuse this matrix because the size didn't match.  making a new one");
					this.rongTransferMatrix = null;
				}
			}

			if (this.rongTransferMatrix == null) { // if we need to make a new one for whatever reason
				this.rongTransferMatrix = this.evaluateTransferMatrix(energyBins, timeBins);
				try {
					CSV.write(rongTransferMatrix, new File("output/transfer matrix.csv"), ','); // save it to disc
				} catch (IOException e) {
					e.printStackTrace();
				}
			}

			if (this.calibrationPrecision == 0) {
				this.trueTransferMatrix = this.rongTransferMatrix;
			}
			else {
				double energyResolutionModifier = 1 + (2*MC_RANDOM.nextDouble() - 1)*calibrationPrecision;
				double timeResolutionModifier = 1 + (2*MC_RANDOM.nextDouble() - 1)*calibrationPrecision;
				System.out.println("augmenting ener1gy resolution by " + energyResolutionModifier + " and time resolution by " + timeResolutionModifier);
				for (int i = 0; i < 4; i++)
					this.cosyMapping.coefficients[i][0] *= energyResolutionModifier;
				for (int i = 0; i < 6; i++)
					if (i != 4)
						this.cosyMapping.coefficients[i][4] *= timeResolutionModifier;
				this.trueTransferMatrix = this.evaluateTransferMatrix(energyBins, timeBins);// evaluateTransferMatrix(); // now make up the actual transfer matrix
				for (int i = 0; i < 4; i++)
					this.cosyMapping.coefficients[i][0] /= energyResolutionModifier;
				for (int i = 0; i < 6; i++)
					if (i != 4)
						this.cosyMapping.coefficients[i][4] /= timeResolutionModifier;
			}
		}
	}

	/**
	 * determine the response of the MRSt to neutrons at particular energies. also tack on some
	 * finite difference laplacian roughness measures.
	 * @param energyBins the neutron energy bins by which we linearize [MeV]
	 * @param timeBins the neutron time bins by which we linearize [ns]
	 */
	private double[][] evaluateTransferMatrix(double[] energyBins, double[] timeBins) {
		int n = (energyBins.length - 1)*(timeBins.length - 1);
		double[][] matrix = new double[n][n];
		int j0 = timeBins.length/2; // time symmetry means we only need to evaluate at one time
		for (int i = 0; i < energyBins.length - 1; i++) { // sweep through all energies
			double energy0 = energyBins[i], energy1 = energyBins[i + 1]; // [MeV]
			double time0 = timeBins[j0]*ns, time1 = timeBins[j0 + 1]*ns; // [s]
			double weight = this.efficiency((energy0 + energy1)/2)/TRANSFER_MATRIX_TRIES;
			for (int k = 0; k < TRANSFER_MATRIX_TRIES; k++) {
				double energyI = energy0 + MC_RANDOM.nextDouble()*(energy1 - energy0); // randomly choose values from the bin [MeV]
				double timeI = time0 + MC_RANDOM.nextDouble()*(time1 - time0); // [s]

				double[] etUncorrected = simulate(energyI, timeI, true);
				if (!Double.isNaN(etUncorrected[0])) { // sometimes, they won't hit the CsI cathode. That's fine.
					double[] et = backCalculate(etUncorrected); // do the simulation!

					double energyO = et[0], timeO = et[1]/ns; // then convert to the same units as the bins
					int eBin = Math2.bin(energyO, energyBins);
					int tBin = Math2.bin(timeO, timeBins);
					if (eBin >= 0 && eBin < energyBins.length - 1 && tBin >= 0 && tBin < timeBins.length - 1) // if it falls in detectable bounds
						matrix[(timeBins.length - 1)*eBin + tBin][(timeBins.length - 1)*i + j0] += weight; // add it to the row
				}
			}
		}

		for (int i = 0; i < n; i++) { // now iterate through all of the rows
			for (int j = 0; j < n; j++) { // and all of the columns, of which we still don't have many
				int jRef = j/(timeBins.length - 1)*(timeBins.length - 1) + j0; // look for the nearest column that is filled
				int iRef = i - j + jRef; // [i,j] is equivalent to an [iRef,jRef] by time symmetry
				if (iRef >= 0 && i/(timeBins.length - 1) == iRef/(timeBins.length - 1)) // this may have jumped over to the next energy,
					matrix[i][j] = matrix[iRef][jRef];
				else
					matrix[i][j] = 0; // but if so, it's almost certainly 0
			}
		}

		return matrix;
	}

	/**
	 * estimate the original time and energy of this ion's neutron without looking at its actual
	 * time and energy, by guessing its energy and accounting for travel time.
	 * @param coords the position and time where it hits the focal plane (m, m, m, s)
	 * @return { energy, time } [MeV, s].
	 */
	private double[] backCalculate(double... coords) {
		double E = energyVsPosition.evaluate(coords[x]/Math.cos(focalPlaneAngle)); // [J]
		double v = Math.sqrt(2*E/cosyMapping.ion.mass);
		double focusingDistance = coords[z];
		double Δt = cosyMapping.getT(0, 0, 0, 0, 0, E/MeV); // re-use the COSY mapping to estimate the time of flight from energy
		return new double[] {
			  E/energyFactor/MeV,
			  coords[3] - (Δt + focusingDistance/v) };
	}

	/**
	 * choose a random location in the foil for the neutron to collide.
	 * @return { x, y, z } [m]
	 */
	private double[] chooseCollisionPosition() {
		double xF = foilWidth/2*(2*MC_RANDOM.nextDouble()-1);
		double yF = foilHeight/2*(2*MC_RANDOM.nextDouble()-1);
		double zF = foilDistance + foilThickness/2*(2*MC_RANDOM.nextDouble()-1); // assume foil is thin, so every z coordinate is equally likely
		return new double[] { xF, yF, zF };
	}

	/**
	 * choose a random location in the aperture plane for the deuteron to pass through.
	 * @return { x, y, z } [m]
	 */
	private double[] chooseAperturePosition() {
		double xA = (2*MC_RANDOM.nextDouble()-1)*apertureWidth/2; // assume aperture is far away, so every point in it is equally likely to be hit
		double yA = (2*MC_RANDOM.nextDouble()-1)*apertureHeight/2;
		double zA = apertureDistance;
		return new double[] { xA, yA, zA };
	}

	/**
	 * compute the velocity with which the deuteron passes through the aperture.
	 * @param energy the energy of the initial particle [MeV]
	 * @param rFoil {x,y,z} of the point at which the neutron strikes the deuteron [m].
	 * @param rAperture {x,y,z} of the point where the deuteron passes through the aperture [m].
	 * @return { vx, vy, vz } [m/s]
	 */
	private double[] computeFinalVelocity(
			double energy, double[] rFoil, double[] rAperture, boolean foilBlur) {
		double E0 = energy*MeV; // convert energy from [MeV] to [J]

		double[] nHat = { // get the unit vector in the direction of the neutron
				rFoil[x], rFoil[y], rFoil[z] };
		double norm = Math.sqrt(Math2.sqr(nHat));
		for (int i = 0; i < 3; i ++)
			nHat[i] /= norm;

		double[] dHat = { // and the unit vector in the direction of the ion
				rAperture[x] - rFoil[x], rAperture[y] - rFoil[y], rAperture[z] - rFoil[z] };
		norm = Math.sqrt(Math2.sqr(dHat));
		for (int i = 0; i < 3; i ++)
			dHat[i] /= norm;

		double E1;
		if (foilBlur) {
			double cosθ = nHat[x]*dHat[x] + nHat[y]*dHat[y] + nHat[z]*dHat[z];
			E1 = energyFactor*E0*cosθ*cosθ; // assume elastic collision between neutron and ion
			double distance = (foilDistance + foilThickness/2 - rFoil[z])/dHat[z];
			E1 = energyVsDistance.evaluate(distanceVsEnergy.evaluate(E1) - distance); // lose some energy to stopping in the foil
		}
		else {
			E1 = energyFactor*E0;
		}

		if (E1 < minEd || E1 > maxEd) {
			return new double[] { Double.NaN, Double.NaN, Double.NaN }; // some won't make it through the "energy aperture"
		}
		else {
			double v = Math.sqrt(2*E1/cosyMapping.ion.mass); // get the resulting velocity otherwise
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

		double KInit = 1/2.*cosyMapping.ion.mass*Math2.sqr(vInit)/MeV; // for the 'd' coordinate, we must convert velocity to energy
		double xPlane = cosyMapping.getX(rFoil[x], vInit[x], rFoil[y], vInit[y], tNeutron, KInit);
		double vxPlane = cosyMapping.getVx(rFoil[x], vInit[x], rFoil[y], vInit[y], tNeutron, KInit);
		double yPlane = cosyMapping.getY(rFoil[x], vInit[x], rFoil[y], vInit[y], tNeutron, KInit);
		double vyPlane = cosyMapping.getVy(rFoil[x], vInit[x], rFoil[y], vInit[y], tNeutron, KInit);
		double[] rPlane = { xPlane, yPlane, 0 }; // [m]
		double[] vFinal = { vxPlane, vyPlane, vInit[z] }; // [m/s]
		double tPlane = cosyMapping.getT(rFoil[x], vInit[x], rFoil[y], vInit[y], tNeutron, KInit);

		double tFocusing = rPlane[x] / (vFinal[z]/Math.tan(focalPlaneAngle) - vFinal[x]); // finally, account for the additional t that it takes to strike the focal plane
		double[] rFocused = {
				rPlane[x] + tFocusing*vFinal[x],
				rPlane[y] + tFocusing*vFinal[y],
				rPlane[z] + tFocusing*vFinal[z] };

		return new double[] { rFocused[x], rFocused[y], rFocused[z], tPlane + tFocusing };
	}
}
