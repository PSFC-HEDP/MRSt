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
	
	private final double probHitsFoil; // probability that the neutron hits the foil
	
	/**
	 * perform some preliminary calculations for the provided configuration.
	 */
	public MonteCarlo(double foilDistance, double foilRadius, double foilThickness) {
		this.foilDistance = foilDistance;
		this.foilThickness = foilThickness;
		
		double foilMaxAngle = Math.atan(foilRadius/foilDistance);
		this.probHitsFoil = (1 - Math.cos(foilMaxAngle))/2;
	}
	
	/**
	 * simulate a single random neutron emitted from TCC at the given energy and determine the position and velocity at which its child deuteron crosses the focal plane.
	 * @param energy initial energy of released neutron [eV].
	 * @return {x, y, z, vx, vy, vz} [m].
	 */
	public double[] response(double energy) {
		System.out.print("[");
		
		double[] rCollision = chooseCollisionPosition();
		System.out.print(String.format("%f, %f, %f,", rCollision[0], rCollision[1], rCollision[2]));
		
		double[] vInitial = computeInitialVelocity(energy, rCollision);
		
		double[] rAperture = chooseAperturePosition();
		
		double[] vFinal = computeFinalVelocity(vInitial, rCollision, rAperture);
		
		double[] rFocal = computeFocusedPosition(rAperture, vFinal);
		
		System.out.println("],");
		return rFocal;
	}
	
	/**
	 * choose a random location in the foil for the neutron to collide.
	 * @return {x, y, z}
	 */
	private double[] chooseCollisionPosition() {
		double θ = Math.acos(1 - 2*Math.random()*probHitsFoil);
		double r = foilDistance*Math.tan(θ); // NOTE: original code assumes uniform distribution within foil; I account for nonzero solid angle subtended at TCC.
		double φ = Math.random()*2*Math.PI;
		double z = foilDistance + (2*Math.random()-1)*foilThickness/2; // assume foil is thin, so every z coordinate is equally likely
		return new double[] { r*Math.cos(φ), r*Math.sin(φ), z };
	}

	private double[] computeInitialVelocity(double energy, double[] rCollision) {
		// TODO: Implement this
		return null;
	}

	private double[] chooseAperturePosition() {
		// TODO: Implement this
		return null;
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
		MonteCarlo sim = new MonteCarlo(3.0e-3, 3.0e-4, 80e-6);
		for (int i = 0; i < 10000; i ++)
			sim.response(14.2e6);
//		System.out.println(Arrays.toString(sim.response(14)));
	}
	
}
