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

import java.util.Locale;

/**
 * a file with some useful numerical analysis stuff.
 * 
 * @author Justin Kunimune
 */
public class NumericalMethods {
	
	/**
	 * draw a number from a Gaussian distribution.
	 * @param μ mean
	 * @param σ standard deviation
	 * @return the number
	 */
	public static double normal(double μ, double σ) {
		if (σ < 0)
			throw new IllegalArgumentException("standard deviation must not be negative");
		double u1 = Math.random();
		double u2 = Math.random();
		double z = Math.sqrt(-2*Math.log(u1))*Math.cos(2*Math.PI*u2);
		return σ*z + μ;
	}
	
	/**
	 * draw a number from a Poisson distribution.
	 * @param λ expectation value
	 * @return the number
	 */
	public static int poisson(double λ) {
		if (λ < 20)
			throw new IllegalArgumentException("this approximation is bad for expectation < 20");
		return (int) Math.round(Math.max(0, normal(λ, Math.sqrt(λ))));
	}
	
	/**
	 * a discrete representation of an unknown function, capable of evaluating in log time.
	 * 
	 * @author Justin Kunimune
	 */
	public static class DiscreteFunction {
		
		private final int resolution; // either the number of equal intervals in the x array, or 0 to indicate unequal intervals
		private double[] X, Y;
		
		/**
		 * instantiate a new function given raw data. x must monotonically increase, or the
		 * evaluation method won't work.
		 * @param x the x values
		 * @param y the corresponding y values
		 */
		public DiscreteFunction(double[] x, double[] y) {
			if (x.length != y.length)
				throw new IllegalArgumentException("datum lengths must match");
			for (int i = 1; i < x.length; i ++)
				if (x[i] < x[i-1])
					throw new IllegalArgumentException("x must be monotonically increasing.");
			this.X = x;
			this.Y = y;
			this.resolution = 0;
		}
		
		/**
		 * instantiate a new function give x data and y data, and assuming x values are all
		 * equally spaced.
		 * @param x the x values had better be properly spaced, because I don't have a good way
		 * to check.
		 * @param y the corresponding y values
		 * @param resolution the number of x intervals
		 */
		public DiscreteFunction(double[] x, double[] y, int resolution) {
			if (x.length != y.length)
				throw new IllegalArgumentException("datums lengths must match");
			for (int i = 1; i < x.length; i ++)
				if (x[i] < x[i-1])
					throw new IllegalArgumentException("x must be monotonically increasing.");
			this.X = x;
			this.Y = y;
			this.resolution = resolution;
		}
		
		/**
		 * it's a function. evaluate it. if this function's x values are equally spaced, this
		 * can be run in O(1) time. otherwise, it will take O(log(n)).
		 * @param x
		 * @return f(x)
		 */
		public double evaluate(double x) {
			int i; // we will linearly interpolate x from (X[i], X[i+1]) onto (Y[i], Y[i+1]).
			if (x < X[0]) // if it's out of bounds, we will extrapolate from the lowest values
				i = 0;
			else if (x > X[X.length-1]) // or highest values, depending on which is appropriate
				i = X.length-2;
			else if (this.resolution > 0) // nonzero resolution means we can find i itself with linear interpolation
				i = (int)((x - X[0])/(X[resolution] - X[0])*resolution); // linearly interpolate x from X to i
			else { // otherwise, we'll need a binary search
				int min = 0, max = X.length; // you know about binary searches, right?
				i = (min + max)/2;
				while (max - min > 1) { // I probably don't need to explain this.
					if (X[i] < x)
						min = i;
					else
						max = i;
					i = (min + max)/2;
				}
			}
			return (x - X[i])/(X[i+1] - X[i])*(Y[i+1] - Y[i]) + Y[i]; // linearly interpolate x from X[i] to Y[i]

		}
		
		/**
		 * return the inverse of this, assuming it has an increasing inverse.
		 * @return the inverse.
		 */
		public DiscreteFunction inv() {
			try {
				return new DiscreteFunction(this.Y, this.X);
			} catch (IllegalArgumentException e) {
				throw new IllegalArgumentException("cannot invert a non-monotonically increasing function.");
			}
		}
		
		/**
		 * return the antiderivative of this, with some arbitrary vertical shift applied.
		 * @return the antiderivative.
		 */
		public DiscreteFunction antiderivative() {
			double[] yOut = new double[X.length];
			yOut[0] = 0; // arbitrarily set the zeroth point to 0
			for (int i = 1; i < X.length; i ++) {
				yOut[i] = yOut[i-1] + (Y[i-1] + Y[i])/2*(X[i] - X[i-1]); // solve for subsequent points using a trapezoid rule
			}
			return new DiscreteFunction(X, yOut); // NOTE: original code constructs a spline and then integrates that with RK45; I use piecewise lines and integrate exactly
		}
		
		/**
		 * return a copy of this that can be evaluated in O(1) time. some information will be
		 * lost depending on resolution.
		 * @param resolution the desired resolution of the new function.
		 * @return the indexed function.
		 */
		public DiscreteFunction indexed(int resolution) {
			double[] xOut = new double[resolution+1];
			double[] yOut = new double[resolution+1];
			for (int i = 0; i <= resolution; i ++) {
				xOut[i] = (double)i/resolution*(X[X.length-1] - X[0]) + X[0]; // first, linearly create the x on which we are to get solutions
				yOut[i] = this.evaluate(xOut[i]); // then get the y values
			}
			
			return new DiscreteFunction(xOut, yOut, resolution);
		}
		
		@Override
		public String toString() {
			String s = "np.array([";
			for (int i = 0; i < X.length; i ++)
				s += String.format(Locale.US, "[%g,%g],", X[i], Y[i]);
			s += "])";
			return s;
		}
		
	}
	
}
