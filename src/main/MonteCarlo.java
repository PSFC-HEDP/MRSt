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

import java.io.File;
import java.io.IOException;
import java.util.Arrays;

import main.NumericalMethods.DiscreteFunction;

/**
 * the class where all the math is.
 * 
 * @author Justin Kunimune
 */
public class MonteCarlo {
	
	private static final int STOPPING_DISTANCE_RESOLUTION = 36;
	
	private final double foilDistance; // z coordinate of midplane of foil [m]
	private final double foilThickness; // thickness of foil [m]
	private final double apertureDistance; // distance from TCC to aperture [m]
	private final double apertureWidth; // horizontal dimension of aperture [m]
	private final double apertureHeight; // vertical dimension of aperture [m]
	
	private final double probHitsFoil; // probability that the neutron goes through the foil
	private final double probMakesDeuteron; // probability that the neutron interacts with the foil and releases a deuteron
	private final double probHitsAperture; // probability that a secondary deuteron goes through the aperture
	
	private final DiscreteFunction distanceVsEnergy;
	private final DiscreteFunction energyVsDistance;
	
	/**
	 * perform some preliminary calculations for the provided configuration.
	 * @param stoppingPower a double[][] containing two columns and n rows. the zeroth column is
	 * the reference values of E in [keV] and the last column is the corresponding values of
	 * dE/dx in [keV/μm].
	 */
	public MonteCarlo(
			double foilDistance, double foilRadius, double foilThickness,
			double foilDensity, double foilCrossSection, double[][] stoppingPowerData,
			double apertureDistance, double apertureWidth, double apertureHeight) {
		this.foilDistance = foilDistance;
		this.foilThickness = foilThickness;
		this.apertureDistance = apertureDistance;
		this.apertureWidth = apertureWidth;
		this.apertureHeight = apertureHeight;
		
		double foilMaxAngle = Math.atan(foilRadius/foilDistance);
		this.probHitsFoil = (1 - Math.cos(foilMaxAngle))/2;
		this.probMakesDeuteron = foilDensity*foilCrossSection*foilThickness; // TODO: this should be a function of energy
		this.probHitsAperture = apertureWidth*apertureHeight /
				(4*Math.PI*Math.pow(apertureDistance,2)); // TODO: account for anisotropic scattering
		
		double[] dxdE = new double[stoppingPowerData.length]; // integrate the stopping power to get stopping distance
		double[] E = new double[stoppingPowerData.length];
		for (int i = 0; i < stoppingPowerData.length; i ++) {
			dxdE[i] = 1/(stoppingPowerData[i][1]*1e9*(-Particle.E.charge)); // converting from [GeV/m] to [m/J]
			E[i] = stoppingPowerData[i][0]*1e3*(-Particle.E.charge); // and from [keV] to [J]
		}
		DiscreteFunction distanceVsEnergyRaw = new DiscreteFunction(E, dxdE).antiderivative();
		this.distanceVsEnergy = distanceVsEnergyRaw.indexed(STOPPING_DISTANCE_RESOLUTION);
		this.energyVsDistance = distanceVsEnergyRaw.inv().indexed(STOPPING_DISTANCE_RESOLUTION);
	}
	
	public double efficiency(double energy) {
		return probHitsFoil * probMakesDeuteron * probHitsAperture;
	}
	
	/**
	 * simulate a single random neutron emitted from TCC at the given energy and determine the position and velocity at which its child deuteron crosses the focal plane.
	 * @param energy initial energy of released neutron [eV].
	 * @param ion either Particle.P or Particle.D depending on the type of foil.
	 * @return {x, y, z, vx, vy, vz} [m].
	 */
	public double[] response(double energy, Particle ion) {
		System.out.print("[");
		
		double[] rCollision = chooseCollisionPosition();
		
		double[] rAperture = chooseAperturePosition();
		
		double[] vFinal = computeFinalVelocity(energy, ion, rCollision, rAperture);
		double v2 = vFinal[0]*vFinal[0] + vFinal[1]*vFinal[1] + vFinal[2]*vFinal[2];
		double E = 1/2.*ion.mass*v2;
		System.out.print(-E/Particle.E.charge);
		
		double[] rFocal = computeFocusedPosition(rAperture, vFinal);
		
		System.out.println("],");
		return rFocal;
	}
	
	/**
	 * choose a random location in the foil for the neutron to collide.
	 * @return { x, y, z } [m]
	 */
	private double[] chooseCollisionPosition() {
		double θ = Math.acos(1 - 2*Math.random()*probHitsFoil);
		double r = foilDistance*Math.tan(θ); // NOTE: original code assumes uniform distribution within foil; I account for nonzero solid angle subtended at TCC.
		double φ = Math.random()*2*Math.PI;
		double z = foilDistance + (2*Math.random()-1)*foilThickness/2; // assume foil is thin, so every z coordinate is equally likely
		return new double[] { r*Math.cos(φ), r*Math.sin(φ), z };
	}
	
	/**
	 * choose a random location in the aperture plane for the deuteron to pass through.
	 * @return { x, y, z } [m]
	 */
	private double[] chooseAperturePosition() {
		double x = (2*Math.random()-1)*apertureWidth/2; // assume aperture is far away, so every point in it is equally likely to be hit
		double y = (2*Math.random()-1)*apertureHeight/2;
		double z = apertureDistance;
		return new double[] { x, y, z };
	}
	
	/**
	 * compute the velocity with which the deuteron passes through the aperture.
	 * @param vInitial {vx,vy,vz} of the neutron as it enters the foil [m/s].
	 * @param A ratio of charged particle mass to neutron mass.
	 * @param rFoil {x,y,z} of the point at which the neutron strikes the deuteron [m].
	 * @param rAperture {x,y,z} of the point where the deuteron passes through the aperture [m].
	 * @return { vx, vy, vz } [m]
	 */
	private double[]
			computeFinalVelocity(double energy, Particle ion, double[] rFoil, double[] rAperture) {
		double E0 = -Particle.E.charge*energy; // convert energy from [eV] to [J]
		double d0x = rFoil[0];
		double d0y = rFoil[1];
		double d0z = rFoil[2];
		double norm = Math.sqrt(d0x*d0x + d0y*d0y + d0z*d0z);
		d0x /= norm;
		d0y /= norm;
		d0z /= norm;
		
		double d1x = rAperture[0] - rFoil[0];
		double d1y = rAperture[1] - rFoil[1];
		double d1z = rAperture[2] - rFoil[2];
		norm = Math.sqrt(d1x*d1x + d1y*d1y + d1z*d1z);
		d1x /= norm;
		d1y /= norm;
		d1z /= norm;
		
		double cosθ = d0x*d1x + d0y*d1y + d0z*d1z;
		double A = ion.mass/Particle.N.mass; // assume elastic collision between neutron and ion
		double E1 = 4*A/Math.pow(A + 1, 2)*E0*cosθ*cosθ; // TODO this is a little redundant; see if it's worth computing ahead of time
		double distance = (foilDistance + foilThickness/2 - rFoil[2])/d0z;
		E1 = energyVsDistance.evaluate(distanceVsEnergy.evaluate(E1) - distance); // lose some energy by dragging through the foil
		double v = Math.sqrt(2*E1/ion.mass); // get the resulting velocity
		return new double[] { v*d1x, v*d1y, v*d1z };
	}
	
	private double[] computeFocusedPosition(double[] rAperture, double[] vFinal) {
		// TODO: Implement this
		return null;
	}
	
	/**
	 * @param args
	 * @throws IOException 
	 * @throws NumberFormatException 
	 */
	public static void main(String[] args) throws NumberFormatException, IOException {
		MonteCarlo sim = new MonteCarlo(
				3.0e-3, 0.3e-3, 80e-6,
				10, 10, CSV.read(new File("data/stopping_power_deuterons.csv"), ','),
				6e0, 4.0e-3, 20.0e-3);
		for (int i = 0; i < 10; i ++)
			sim.response(14e6, Particle.D);
//		System.out.println(Arrays.toString(sim.response(14)));
	}
	
}
