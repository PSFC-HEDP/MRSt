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
package util;

import physics.Particle;

/**
 * a table of polynomials output by COSY INFINITY to quasilinearly map one
 * vector space to another
 */
public class COSYMapping {
	private static final double SPEED_OF_LIGHT = 2.99792458e8;

	public final double[][] coefficients;
	public final int[][] exponents;
	public Particle ion;
	public double K0; // central energy (J)
	public double v0; // central velocity (m/s)
	public double t0; // central time (s)
	public double t1; // time scaling thing (s)

	/**
	 * create a new COSY map
	 * @param coefficients the coefficient table; each colum is a polynomial
	 * @param exponents the exponent table; each colum is an input variable
	 */
	public COSYMapping(double[][] coefficients, int[][] exponents, Particle ion, double neutronEnergy) {
		this.coefficients = coefficients;
		this.exponents = exponents;
		this.ion = ion;
		double A = ion.mass/Particle.N.mass;
		this.K0 = neutronEnergy*4*A/Math.pow(1 + A, 2);
		double K0_joule = K0*1e6*(-Particle.E.charge);
		this.v0 = Math.sqrt(2*K0_joule/ion.mass); // and get the corresponding speed
		double γ = Math.pow(1 - Math.pow(v0/SPEED_OF_LIGHT, 2), -1/2.);
		double L = 8; // this number doesn't matter so I'm eyeballing it
		this.t0 = L/v0; // and corresponding time
		this.t1 = -(1+γ)/γ/v0; // and why is time measured in meters?
	}

	private double[] cleanInput(
		  double x, double vx, double y, double vy, double t, double K) {
		return new double[] {x, vx/v0, y, vy/v0, (t - t0)/t1, (K - K0)/K0}; // COSY takes velocities in [rad] and times in [m]
	}

	public double getX(double x, double vx, double y, double vy, double t, double K) {
		return this.polynomial(0, cleanInput(x, vx, y, vy, t, K));
	}

	public double getVx(double x, double vx, double y, double vy, double t, double K) {
		return this.polynomial(1, cleanInput(x, vx, y, vy, t, K))*v0;
	}

	public double getY(double x, double vx, double y, double vy, double t, double K) {
		return this.polynomial(2, cleanInput(x, vx, y, vy, t, K));
	}

	public double getVy(double x, double vx, double y, double vy, double t, double K) {
		return this.polynomial(3, cleanInput(x, vx, y, vy, t, K))*v0;
	}

	public double getT(double x, double vx, double y, double vy, double t, double K) {
		return this.polynomial(4, cleanInput(x, vx, y, vy, t, K))*t1 + t0;
	}

	/**
	 * evaluate a single polynomial from that table of coefficients.
	 * @param i the index of the parameter to compute
	 * @param input the initial values to use
	 * @return the final value of parameter i
	 */
	private double polynomial(int i, double... input) {
		double output = 0;
		for (int j = 0; j < coefficients.length; j ++) {
			double term = coefficients[j][i];
			for (int k = 0; k < input.length; k ++) {
				term *= Math.pow(input[k], exponents[j][k]);
			}
			output += term;
		}
		return output;
	}

}
