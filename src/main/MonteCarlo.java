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
	
	/**
	 * perform some preliminary calculations for the provided configuration.
	 */
	public MonteCarlo() {
		
	}
	
	/**
	 * simulate a single random neutron emitted from TCC at the given energy and determine the position and velocity at which it crosses the focal plane.
	 * @param energy in MeV.
	 * @return {x, y, z, vx, vy, vz} in m.
	 */
	public double[] response(double energy) {
		return new double[] {0, 0, 0, 0, 0, 0};
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		MonteCarlo sim = new MonteCarlo();
		System.out.println(Arrays.toString(sim.response(14)));
	}
	
}
