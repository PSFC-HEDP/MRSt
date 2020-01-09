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

import java.util.Arrays;

/**
 * the class where all the math is.
 * 
 * @author Justin Kunimune
 */
public class MonteCarlo {
	
	private final double foilDistance; // z coordinate of midplane of foil [m]
	private final double foilThickness; // thickness of foil [m]
	private final double apertureDistance; // distance from TCC to aperture [m]
	private final double apertureWidth; // horizontal dimension of aperture [m]
	private final double apertureHeight; // vertical dimension of aperture [m]
	
	private final double probHitsFoil; // probability that the neutron goes through the foil
	private final double probMakesDeuteron; // probability that the neutron interacts with the foil and releases a deuteron
	private final double probHitsAperture; // probability that a secondary deuteron goes through the aperture
	
	/**
	 * perform some preliminary calculations for the provided configuration.
	 */
	public MonteCarlo(
			double foilDistance, double foilRadius, double foilThickness,
			double foilDensity, double foilCrossSection,
			double apertureDistance, double apertureWidth, double apertureHeight) {
		this.foilDistance = foilDistance;
		this.foilThickness = foilThickness;
		this.apertureDistance = apertureDistance;
		this.apertureWidth = apertureWidth;
		this.apertureHeight = apertureHeight;
		
		double foilMaxAngle = Math.atan(foilRadius/foilDistance);
		this.probHitsFoil = (1 - Math.cos(foilMaxAngle))/2;
		this.probMakesDeuteron = foilDensity*foilCrossSection*foilThickness; // TODO: this should be a function of energy
		this.probHitsAperture = apertureWidth*apertureHeight / (4*Math.PI*Math.pow(apertureDistance,2));
	}
	
	public double efficiency(double energy) {
		return probHitsFoil * probMakesDeuteron * probHitsAperture;
	}
	
	/**
	 * simulate a single random neutron emitted from TCC at the given energy and determine the position and velocity at which its child deuteron crosses the focal plane.
	 * @param energy initial energy of released neutron [eV].
	 * @return {x, y, z, vx, vy, vz} [m].
	 */
	public double[] response(double energy) {
		System.out.print("[");
		
		double[] rCollision = chooseCollisionPosition();
		System.out.print(String.format("%f, %f, %f, ", rCollision[0], rCollision[1], rCollision[2]));
		
		double[] vInitial = computeInitialVelocity(energy, rCollision);
		System.out.print(String.format("%f, %f, %f, ", vInitial[0], vInitial[1], vInitial[2]));

		double[] rAperture = chooseAperturePosition();
		System.out.print(String.format("%f, %f, %f, ", rAperture[0], rAperture[1], rAperture[2]));
		
		double[] vFinal = computeFinalVelocity(vInitial, rCollision, rAperture);
		
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
	 * compute the initial velocity of a neutron given the point at which it hits the foil and its initial energy.
	 * @param energy initial energy of neutron [eV].
	 * @param rCollision {x,y,z} at which it dies in the foil [m].
	 * @return { vx, vy, vz } [m/s]
	 */
	private double[] computeInitialVelocity(double energy, double[] rCollision) {
		energy = -Particle.E.charge*energy; // convert energy from eV to J
		double vx = rCollision[0]; // assume emitted from TCC
		double vy = rCollision[1];
		double vz = rCollision[2];
		double r = Math.sqrt(vx*vx + vy*vy + vz*vz);
		double v = Math.sqrt(2*energy/Particle.N.mass); // assume nonrelativistic
		vx *= v/r;
		vy *= v/r;
		vz *= v/r;
		return new double[] { vx, vy, vz };
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
	
	private double[]
			computeFinalVelocity(double[] vInitial, double[] rCollision, double[] rAperture) {
		// TODO: Implement this
		return null;
	}

	private double[] computeFocusedPosition(double[] rAperture, double[] vFinal) {
		// TODO: Implement this
		return null;
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		MonteCarlo sim = new MonteCarlo(
				3.0e-3, 0.3e-3, 80e-6,
				10, 10,
				6e0, 4.0e-3, 20.0e-3);
		for (int i = 0; i < 1000; i ++)
			sim.response(14.2e6);
//		System.out.println(Arrays.toString(sim.response(14)));
	}
	
}
