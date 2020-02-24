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
		if (λ < 20) {
			double u = Math.random()*Math.exp(λ);
			long kFact = 1;
			for (int k = 0; k < 40; k ++) {
				if (k != 0) kFact *= k;
				u -= Math.pow(λ, k)/kFact;
				if (u < 0)
					return k;
			}
			return 40;
		}
		else { // use Gaussian approximation for high expectations
			return (int) Math.round(Math.max(0, normal(λ, Math.sqrt(λ))));
		}
	}
	
	/**
	 * compute the nth moment of the histogram over the whole domain. normalize and center it,
	 * if applicable.
	 * @param x the bin edges
	 * @param y the number in each bin
	 * @return the nth [normalized] [centered] [normalized] moment
	 */
	public static double moment(int n, double[] x, double[] y) {
		return moment(n, x, y, x[0], x[x.length-1]);
	}
	
	/**
	 * compute the nth moment of the histogram. normalize and center it, if applicable.
	 * @param x the bin edges
	 * @param y the number in each bin
	 * @param a the lower integration bound
	 * @param b the upper integration bound
	 * @return the nth [normalized] [centered] [normalized] moment
	 */
	public static double moment(int n, double[] x, double[] y, double a, double b) {
		if (x.length != y.length+1)
			throw new IllegalArgumentException("Array lengths do not correspond.");
		double N = (n > 0) ? moment(0, x, y, a, b) : 1;
		double μ = (n > 1) ? moment(1, x, y, a, b) : 0;
		double σ = (n > 2) ? Math.sqrt(moment(2, x, y, a, b)) : 1;
		double sum = 0;
		for (int i = 0; i < y.length; i ++) {
			double xL = Math.max(a, x[i]); // define the bounds of this integrand bin, which might not be the bounds of the datum bin
			double xR = Math.min(b, x[i+1]);
			double w = (xR - xL)/(x[i+1] - x[i]); // determine what fraction of the data in this bin fall into the integrand bin
			sum += w*y[i]*Math.pow(((xL + xR)/2 - μ)/σ, n); // sum up the average x value to whatever power
		}
		return sum/N;
	}
	
	/**
	 * compute the 0th moment of the histogram
	 * @param x the bin edges
	 * @param y the number in each bin
	 * @return the total number of things counted
	 */
	public static double definiteIntegral(double[] x, double[] y) {
		return definiteIntegral(x, y, x[0], x[x.length-1]);
	}
	
	/**
	 * compute the 0th moment of the histogram
	 * @param x the bin edges
	 * @param y the number in each bin
	 * @param a the lower integration bound
	 * @param b the upper integration bound
	 * @return the total number of things counted
	 */
	public static double definiteIntegral(double[] x, double[] y, double a, double b) {
		if (x.length != y.length+1)
			throw new IllegalArgumentException("Array lengths do not correspond.");
		double s = 0;
		for (int i = 0; i < y.length; i ++) {
			double wl = Math.max(0, Math.min(1, (x[i+1] - a)/(x[i+1] - x[i])));
			double wr = Math.max(0, Math.min(1, (b - x[i])/(x[i+1] - x[i])));
			s += (wl+wr-1)*y[i];
		}
		return s;
	}
	
	/**
	 * compute the mean of the histogram
	 * @param x the bin edges
	 * @param y the number in each bin
	 * @return the normalized 1st moment
	 */
	public static double mean(double[] x, double[] y) {
		return moment(1, x, y);
	}
	
	/**
	 * compute the mean of the histogram
	 * @param x the bin edges
	 * @param y the number in each bin
	 * @param a the lower integration bound
	 * @param b the upper integration bound
	 * @return the normalized 1st moment
	 */
	public static double mean(double[] x, double[] y, double a, double b) {
		return moment(1, x, y, a, b);
	}
	
	/**
	 * compute the standard deviation of the histogram
	 * @param x the bin edges
	 * @param y the number in each bin
	 * @return the square root of the normalized centered 2nd moment
	 */
	public static double std(double[] x, double[] y) {
		return Math.sqrt(moment(2, x, y));
	}
	
	/**
	 * compute the standard deviation of the histogram
	 * @param x the bin edges
	 * @param y the number in each bin
	 * @param a the lower integration bound
	 * @param b the upper integration bound
	 * @return the square root of the normalized centered 2nd moment
	 */
	public static double std(double[] x, double[] y, double a, double b) {
		return Math.sqrt(moment(2, x, y, a, b));
	}
	
	/**
	 * find the last index of the highest value
	 * @param x the array of values
	 * @return i such that x[i] >= x[j] for all j
	 */
	public static int argmax(double[] x) {
		int argmax = -1;
		for (int i = 0; i < x.length; i ++)
			if (!Double.isNaN(x[i]) && (argmax == -1 ||
					x[i] > x[argmax]))
				argmax = i;
		return argmax;
	}
	
	/**
	 * find the second order finite difference derivative. for best results, x
	 * should be evenly spaced.
	 * @param x the x values
	 * @param y the corresponding y values
	 * @return the slope dy/dx at each point
	 */
	public static double[] derivative(double[] x, double[] y) {
		if (x.length != y.length)
			throw new IllegalArgumentException("Array lengths do not correspond.");
		double[] dydx = new double[x.length];
		dydx[0] = (-1.5*y[0] + 2.0*y[1] - 0.5*y[2]) /
				(x[2] - x[1]);
		for (int i = 1; i < x.length-1; i ++) {
			dydx[i] = (y[i+1] - y[i-1]) / (x[i+1] - x[i-1]);
		}
		dydx[x.length-1] = (-1.5*y[x.length-1] + 2.0*y[x.length-2] - 0.5*y[x.length-3]) /
				(x[x.length-3] - x[x.length-2]);
		return dydx;
	}
	
	/**
	 * coerce x into the range [min, max]
	 * @param min inclusive minimum
	 * @param max inclusive maximum
	 * @param x floating point value
	 * @return int in the range [min, max]
	 */
	public static int coerce(int min, int max, double x) {
		if (x <= min)
			return min;
		else if (x >= max)
			return max;
		else
			return (int) x;
	}
	
	/**
	 * multiply a vector by a matrix
	 * @param A matrix
	 * @param u vector
	 * @return A.u vector
	 */
	public static double[] matmul(double[][] A, double[] v) {
		if (A[0].length != v.length)
			throw new IllegalArgumentException("Multiply a "+A[0].length+"×"+A.length+" by a "+v.length+"×1? Don't you know how matrix multiplication works?");
		double[] u = new double[A.length];
		for (int i = 0; i < A.length; i ++)
			for (int j = 0; j < A[i].length; j ++)
				u[i] += A[i][j]*v[j];
		return u;
	}
	
	/**
	 * find the inverse of a matrix
	 * @param A matrix
	 * @return A^-1 matrix
	 */
	public static double[][] matinv(double a[][]) {
		int n = a.length;
		double x[][] = new double[n][n];
		double b[][] = new double[n][n];
		int index[] = new int[n];
		for (int i = 0; i < n; i ++)
			b[i][i] = 1;

		gaussian(a, index); // Transform the matrix into an uppre triangle

		for (int i = 0; i < n - 1; ++i) // Update the matrix b[i][j] with the ratios stored
			for (int j = i + 1; j < n; ++j)
				for (int k = 0; k < n; ++k)
					b[index[j]][k] -= a[index[j]][i] * b[index[i]][k];

		for (int i = 0; i < n; ++i) { // Perform backward substitutions
			x[n - 1][i] = b[index[n - 1]][i] / a[index[n - 1]][n - 1];
			for (int j = n - 2; j >= 0; j --) {
				x[j][i] = b[index[j]][i];
				for (int k = j + 1; k < n; k ++) {
					x[j][i] -= a[index[j]][k] * x[k][i];
				}
				x[j][i] /= a[index[j]][j];
			}
		}
		return x;
	}
	
	/**
	 * method to carry out the partial-pivoting Gaussian elimination.
	 * @param a
	 * @param index stores pivoting order
	 */
	public static void gaussian(double a[][], int index[]) {
		int n = index.length;
		double c[] = new double[n];

		// Initialize the index
		for (int i = 0; i < n; ++i)
			index[i] = i;

		// Find the rescaling factors, one from each row
		for (int i = 0; i < n; ++i) {
			double c1 = 0;
			for (int j = 0; j < n; ++j) {
				double c0 = Math.abs(a[i][j]);
				if (c0 > c1)
					c1 = c0;
			}
			c[i] = c1;
		}

		// Search the pivoting element from each column
		int k = 0;
		for (int j = 0; j < n - 1; ++j) {
			double pi1 = 0;
			for (int i = j; i < n; ++i) {
				double pi0 = Math.abs(a[index[i]][j]);
				pi0 /= c[index[i]];
				if (pi0 > pi1) {
					pi1 = pi0;
					k = i;
				}
			}

			// Interchange rows according to the pivoting ordre
			int itmp = index[j];
			index[j] = index[k];
			index[k] = itmp;
			for (int i = j + 1; i < n; ++i) {
				double pj = a[index[i]][j] / a[index[j]][j];

				// Record pivoting ratios below the diagonal
				a[index[i]][j] = pj;

				// Modify othre elements accordingly
				for (int l = j + 1; l < n; ++l)
					a[index[i]][l] -= pj * a[index[j]][l];
			}
		}
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
