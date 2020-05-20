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
	
	public static double[] unimode(double[] x, double[] params) {
		if (params.length != 7)
			throw new IllegalArgumentException("Number of params must be 7");
		return unimode(x, params[0], params[1], params[2], params[3], params[4], params[5], params[6]);
	}
	
	/**
	 * generate a generic unimodal distribution with a handful of parameters with which to play.
	 * it will be the linear combination of an erf stepping from yL to yR with a generalized
	 * skew normal distribution with height ~yPeak, and skew and kurtosis-adjustment given by
	 * tilt and fatten
	 * @param x the x axis on which to generate it
	 * @param yL the limit for negative x
	 * @param yR the limit for positive x
	 * @param yPeak the signed magnitude of the peak in the middle
	 * @param xPeak the location of the center of the function
	 * @param std the standard deviation of the peak
	 * @param tilt the amount to skew the peak (0 is no skew)
	 * @param fatten the exponent on the peak (2 is normal gaussian)
	 * @return y values corresponding to the x values
	 */
	public static double[] unimode(double[] x, double yL, double yR, double yPeak,
			double xPeak, double std, double tilt, double fatten) {
		final double y0 = yR, yS = (yL - yR), yG = yPeak - (yL + yR)/2;
		double[] y = new double[x.length];
		for (int i = 0; i < x.length; i ++) {
			double ξ = (x[i] - xPeak)/std;
			y[i] = y0 + yS*erfc(ξ)/2 +
					yG*Math.exp(-Math.pow(Math.abs(ξ), 4/fatten)/2)*erfc(-tilt*ξ);
		}
		return y;
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
	 * do a standard deviation of a not histogram. it's just a list of numbers.
	 * @param x
	 * @param y
	 * @param a
	 * @param b
	 * @return
	 */
	public static double std(double[] x) {
		double mean = 0;
		double meanSqr = 0;
		for (int i = 0; i < x.length; i ++) {
			mean += x[i]/x.length;
			meanSqr += Math.pow(x[i], 2)/x.length;
		}
		return Math.sqrt(meanSqr - Math.pow(mean, 2));
	}
	
	public static double sum(double[] arr) {
		double s = 0;
		for (double x: arr)
			s += x;
		return s;
	}
	
	public static double max(double[] arr) {
		double max = Double.NEGATIVE_INFINITY;
		for (double x: arr)
			if (x > max)
				max = x;
		return max;
	}
	
	public static double max(double[][] arr) {
		double max = Double.NEGATIVE_INFINITY;
		for (double[] row: arr)
			for (double x: row)
				if (x > max)
					max = x;
		return max;
	}
	
	/**
	 * find the last index of the highest value
	 * @param x the array of values
	 * @return i such that x[i] >= x[j] for all j
	 */
	public static int argmax(double[] x) {
		int argmax = -1;
		for (int i = 0; i < x.length; i ++)
			if (!Double.isNaN(x[i]) && (argmax == -1 || x[i] > x[argmax]))
				argmax = i;
		return argmax;
	}
	
	/**
	 * find the last index of the second highest value
	 * @param x the array of values
	 * @return i such that x[i] >= x[j] for all j
	 */
	public static int argpenmax(double[] x) {
		int argmax = argmax(x);
		int argpenmax = -1;
		for (int i = 0; i < x.length; i ++)
			if (i != argmax && !Double.isNaN(x[i]) && (argpenmax == -1 || x[i] > x[argpenmax]))
				argpenmax = i;
		return argpenmax;
	}
	
	/**
	 * find the last index of the lowest value
	 * @param x the array of values
	 * @return i such that x[i] >= x[j] for all j
	 */
	public static int argmin(double[] x) {
		int argmin = -1;
		for (int i = 0; i < x.length; i ++)
			if (!Double.isNaN(x[i]) && (argmin == -1 ||
					x[i] < x[argmin]))
				argmin = i;
		return argmin;
	}
	/**
	 * find the interpolative index of the highest value
	 * @param x the array of values
	 * @return i such that x[i] >= x[j] for all j
	 */
	public static double quadargmax(double[] x) {
		int i = argmax(x);
		if (i <= 0 || Double.isNaN(x[i-1]) || i == x.length-1 || Double.isNaN(x[i+1])) return i;
		double dxdi = (x[i+1] - x[i-1])/2;
		double d2xdi2 = (x[i+1] - 2*x[i] + x[i-1]);
		assert d2xdi2 < 0;
		return i - dxdi/d2xdi2;
	}
	
	/**
	 * find the x coordinate of the highest value
	 * @param x the horizontal axis
	 * @param y the array of values
	 * @return x such that y(x) >= y(z) for all z
	 */
	public static double quadargmax(double[] x, double[] y) {
		try {
			return interp(x, quadargmax(y));
		} catch (IndexOutOfBoundsException e) { // y is empty or all NaN
			return -1;
		}
	}
	
	/**
	 * take the floating-point index of an array using linear interpolation.
	 * @param x the array of values
	 * @param i the partial index
	 * @return x[i], more or less
	 */
	public static double interp(double[] x, double i) {
		if (i < 0 || i > x.length-1)
			throw new IndexOutOfBoundsException("Even partial indices have limits: "+i);
		if (i == x.length-1)
			return x[(int)i];
		return (1 - i%1)*x[(int)i] + (i%1)*x[(int)i+1];
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
	 * convert this 2d histogram to a lower resolution. the output bins must be uniform.
	 * @param xI the horizontal bin edges of the input histogram
	 * @param yI the vertical bin edges of the input histogram
	 * @param zI the counts of the input histogram
	 * @param xO the horizontal bin edges of the desired histogram. these must be uniform.
	 * @param yO the vertical bin edges of the desired histogram. these must be uniform.
	 * @return zO the counts of the new histogram
	 */
	public static double[][] downsample(double[] xI, double[] yI, double[][] zI,
			double[] xO, double[] yO) {
		if (yI.length-1 != zI.length || xI.length-1 != zI[0].length)
			throw new IllegalArgumentException("Array sizes don't match fix it.");
		
		double[][] zO = new double[yO.length-1][xO.length-1]; // first resize the input spectrum to match the transfer matrix
		for (int iI = 0; iI < yI.length-1; iI ++) {
			for (int jI = 0; jI < xI.length-1; jI ++) { // for each small pixel on the input spectrum
				double iO = (yI[iI] - yO[0])/(yO[1] - yO[0]); // find the big pixel of the scaled spectrum
				double jO = (xI[jI] - xO[0])/(xO[1] - xO[0]); // that contains the upper left corner
				int iOint = (int)Math.floor(iO);
				double iOmod = iO - iOint;
				int jOint = (int)Math.floor(jO);
				double jOmod = jO - jOint;
				double cU = Math.min(1, (1 - iOmod)*(yO[1] - yO[0])/(yI[iI+1] - yI[iI])); // find the fraction of it that is above the next pixel
				double cL = Math.min(1, (1 - jOmod)*(xO[1] - xO[0])/(xI[jI+1] - xI[jI])); // and left of the next pixel
				
				addIfInBounds(zO, iOint,   jOint,   zI[iI][jI]*cU*cL); // now add the contents of this spectrum
				addIfInBounds(zO, iOint,   jOint+1, zI[iI][jI]*cU*(1-cL)); // being careful to distribute them properly
				addIfInBounds(zO, iOint+1, jOint,   zI[iI][jI]*(1-cU)*cL); // (I used this convenience method because otherwise I would have to check all the bounds)
				addIfInBounds(zO, iOint+1, jOint+1, zI[iI][jI]*(1-cU)*(1-cL));
			}
		}
		
		return zO;
	}
	
	/**
	 * interpolate y0 from x0 to x1 using a 3rd order spline
	 * @param x1 the desired interpolation points
	 * @param x0 the locations of the spline points
	 * @param y0 the values at the spline points
	 * @return y1 the values at the interpolation points
	 */
	public static double[] spline(double[] x1, double[] x0, double[] y0) {
		float[] x = new float[x0.length];
		float[] y = new float[y0.length];
		for (int i = 0; i < x0.length; i ++) {
			x[i] = (float) x0[i];
			y[i] = (float) y0[i];
		}
		Spline s = Spline.createSpline(x, y);
		double[] y1 = new double[x1.length];
		for (int i = 0; i < x1.length; i ++)
			y1[i] = s.interpolate((float) x1[i]);
		return y1;
	}
	
	/**
	 * a simple convenience method to avoid excessive if statements
	 * @param arr
	 * @param i
	 * @param j
	 * @param val
	 */
	private static void addIfInBounds(double[][] arr, int i, int j, double val) {
		if (i >= 0 && i < arr.length)
			if (j >= 0 && j < arr[i].length)
				arr[i][j] += val;
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
	
	
	private static final double[] cof = {
		-1.3026537197817094, 6.4196979235649026e-1,
		1.9476473204185836e-2, -9.561514786808631e-3, -9.46595344482036e-4,
		3.66839497852761e-4, 4.2523324806907e-5, -2.0278578112534e-5,
		-1.624290004647e-6, 1.303655835580e-6, 1.5626441722e-8, -8.5238095915e-8,
		6.529054439e-9, 5.059343495e-9, -9.91364156e-10, -2.27365122e-10,
		9.6467911e-11, 2.394038e-12, -6.886027e-12, 8.94487e-13, 3.13092e-13,
		-1.12708e-13, 3.81e-16, 7.106e-15, -1.523e-15, -9.4e-17, 1.21e-16, -2.8e-17
	};
	
	/**
	 * The Gauss error function.
	 */
	public static double erf(double x) {
		if (x >= 0.) {
			return 1.0 - erfccheb(x);
		} else {
			return erfccheb(-x) - 1.0;
		}
	}
	
	/**
	 * the complementary Gauss error function
	 */
	public static double erfc(double x) {
		return 1 - erf(x);
	}
	
	private static double erfccheb(double z) {
		double t, ty, tmp, d = 0., dd = 0.;
		if (z < 0.) {
			throw new IllegalArgumentException("erfccheb requires nonnegative argument");
		}
		t = 2. / (2. + z);
		ty = 4. * t - 2.;
		for (int j = cof.length - 1; j > 0; j--) {
			tmp = d;
			d = ty * d - dd + cof[j];
			dd = tmp;
		}
		return t * Math.exp(-z * z + 0.5 * (cof[0] + ty * d) - dd);
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
	
	
	public static final void main(String[] args) {
		int n = 40;
		double[] x = new double[n+1];
		for (int i = 0; i < n+1; i ++)
			x[i] = -3 + 6./n*i;
		System.out.println(Arrays.toString(x));
		System.out.println(Arrays.toString(unimode(x, .5, .5, -2, 0, 1, 6, 2)));
		System.out.println(poisson(0));
	}
	
}
