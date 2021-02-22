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
import java.util.Random;

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
		return normal(μ, σ, Math.random(), Math.random());
	}
	
	/**
	 * draw a number from a Gaussian distribution.
	 * @param μ mean
	 * @param σ standard deviation
	 * @param random the rng to use
	 * @return the number
	 */
	public static double normal(double μ, double σ, Random random) {
		return normal(μ, σ, random.nextDouble(), random.nextDouble());
	}
	
	/**
	 * draw a number from a Gaussian distribution.
	 * @param μ mean
	 * @param σ standard deviation
	 * @return the number
	 */
	private static double normal(double μ, double σ, double u1, double u2) {
		if (σ < 0)
			throw new IllegalArgumentException("standard deviation must not be negative");
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
			return poisson(λ, Math.random());
		else
			return (int) Math.round(Math.max(0, normal(λ, Math.sqrt(λ))));
	}
	
	/**
	 * draw a number from a Poisson distribution.
	 * @param λ expectation value
	 * @param random the rng to use
	 * @return the number
	 */
	public static int poisson(double λ, Random random) {
		if (λ < 20)
			return poisson(λ, random.nextDouble());
		else
			return (int) Math.round(Math.max(0, normal(λ, Math.sqrt(λ), random)));
	}
	
	/**
	 * draw a number from a Poisson distribution.
	 * @param λ expectation value
	 * @return the number
	 */
	private static int poisson(double λ, double u) {
		if (λ < 20) {
			u *= Math.exp(λ);
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
			throw new IllegalArgumentException("You should use a Gaussian approximation for high expectations, but I can't do that with this poisson(double, double) call");
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
	 * compute the nth moment of the histogram over the whole domain. normalize and center it,
	 * if applicable.
	 * @param x the bin edges
	 * @param y the number in each bin
	 * @return the nth [normalized] [centered] [normalized] moment
	 */
	public static Quantity moment(int n, double[] x, Quantity[] y) {
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
	public static Quantity moment(int n, double[] x, Quantity[] y, double a, double b) {
		if (x.length != y.length+1)
			throw new IllegalArgumentException("Array lengths do not correspond.");
		Quantity N = (n > 0) ? moment(0, x, y, a, b) : new Quantity(1, y[0].getN());
		Quantity μ = (n > 1) ? moment(1, x, y, a, b) : new Quantity(0, y[0].getN());
		Quantity σ = (n > 2) ? moment(2, x, y, a, b).sqrt() : new Quantity(1, y[0].getN());
		Quantity sum = new Quantity(0, y[0].getN());
		for (int i = 0; i < y.length; i ++) {
			double xL = Math.max(a, x[i]); // define the bounds of this integrand bin, which might not be the bounds of the datum bin
			double xR = Math.min(b, x[i+1]);
			double w = (xR - xL)/(x[i+1] - x[i]); // determine what fraction of the data in this bin fall into the integrand bin
			sum = sum.plus(y[i].times(w).times(μ.minus((xL + xR)/2).times(-1).over(σ).pow(n))); // sum up the average x value to whatever power
		}
		return sum.over(N);
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
	 * compute the averaged value
	 * @param y the values
	 * @param f the weights
	 */
	public static Quantity average(Quantity[] y, Quantity[] f) {
		Quantity p0 = new Quantity(0, y[0].getN());
		Quantity p1 = new Quantity(0, y[0].getN());
		for (int i = 0; i < y.length; i ++) {
			p0 = p0.plus(f[i]);
			p1 = p1.plus(f[i].times(y[i]));
		}
		return p1.over(p0);
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
	
	public static double sum(double[][] arr) {
		double s = 0;
		for (double[] row: arr)
			for (double x: row)
				s += x;
		return s;
	}
	
	public static double mean(double[] arr) {
		return sum(arr)/arr.length;
	}
	
	public static double fwhm(double[] x, double[] y) {
		int max = argmax(y);
		double xR = Double.POSITIVE_INFINITY;
		for (int i = max + 1; i < y.length; i ++) {
			if (y[i] < y[max]/2.) {
				double c = (y[max]/2. - y[i])/(y[i-1] - y[i]);
				xR = (c*x[i+1] + x[i] + (1-c)*x[i-1])/2.;
				break;
			}
		}
		double xL = Double.NEGATIVE_INFINITY;
		for (int i = max; i >= 1; i --) {
			if (y[i-1] < y[max]/2.) {
				double c = (y[max]/2. - y[i])/(y[i-1] - y[i]);
				xL = (c*x[i+1] + x[i] + (1-c)*x[i-1])/2.;
				break;
			}
		}
		return xR - xL;
	}
	
	public static double[] minus(double[] x) {
		double[] out = new double[x.length];
		for (int i = 0; i < out.length; i ++)
			out[i] = -x[i];
		return out;
	}
	
	public static Quantity[] minus(Quantity[] x) {
		Quantity[] out = new Quantity[x.length];
		for (int i = 0; i < out.length; i ++)
			out[i] = x[i].times(-1);
		return out;
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
		return argmax(minus(x));
	}
	
	/**
	 * find the interpolative index of the highest value
	 * @param x the array of values
	 * @return i such that x[i] >= x[j] for all j
	 */
	public static Quantity quadargmin(int left, int right, Quantity[] x) {
		return quadargmax(left, right, minus(x));
	}
	
	/**
	 * find the interpolative index of the highest value
	 * @param x the array of values
	 * @return i such that x[i] >= x[j] for all j
	 */
	public static double quadargmax(double[] x) {
		return quadargmax(0, x.length, x);
	}
	
	/**
	 * find the interpolative index of the highest value in [left, right)
	 * @param left the leftmost acceptable index
	 * @param right the leftmost unacceptable index
	 * @param x the array of values
	 * @return i such that x[i] >= x[j] for all j in [left, right)
	 */
	public static double quadargmax(int left, int right, double[] x) {
		int i = -1;
		for (int j = left; j < right; j ++)
			if (!Double.isNaN(x[j]) && (i == -1 || x[j] > x[i]))
				i = j;
		if (i == left || Double.isNaN(x[i-1]) || i == right-1 || Double.isNaN(x[i+1])) return i;
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
	 * find the interpolative index of the highest value
	 * @param x the array of values
	 * @return i such that x[i] >= x[j] for all j
	 */
	public static Quantity quadargmax(Quantity[] x) {
		return quadargmax(0, x.length, x);
	}
	
	/**
	 * find the interpolative index of the highest value
	 * @param x the array of values
	 * @return i such that x[i] >= x[j] for all j
	 */
	public static Quantity quadargmax(int left, int right, Quantity[] x) {
		int i = -1;
		for (int j = left; j < right; j ++)
			if (i == -1 || x[j].value > x[i].value)
				i = j;
		if (i <= left || i >= right-1)
			return new Quantity(i, x[i].getN());
		Quantity dxdi = (x[i+1].minus(x[i-1])).over(2);
		Quantity d2xdi2 = (x[i+1]).plus(x[i].times(-2)).plus(x[i-1]);
		assert d2xdi2.value < 0;
		return dxdi.over(d2xdi2).times(-1).plus(i);
	}
	
	/**
	 * find the x coordinate of the highest value in the bounds [left, right)
	 * @param left the leftmost acceptable index
	 * @param right the leftmost unacceptable index
	 * @param x the horizontal axis
	 * @param y the array of values
	 * @return x such that y(x) >= y(z) for all z in [x[left], x[right])
	 */
	public static double quadargmax(int left, int right, double[] x, double[] y) {
		if (x.length != y.length)
			throw new IllegalArgumentException("These array lengths don't match.");
		try {
			return interp(x, quadargmax(Math.max(0, left), Math.min(x.length, right), y));
		} catch (IndexOutOfBoundsException e) { // y is empty or all NaN
			return -1;
		}
	}
	
	/**
	 * find the x coordinate of the highest value in the bounds [left, right)
	 * @param left the leftmost acceptable index
	 * @param right the leftmost unacceptable index
	 * @param x the horizontal axis
	 * @param y the array of values
	 * @return x such that y(x) >= y(z) for all z in [x[left], x[right])
	 */
	public static Quantity quadargmax(int left, int right, double[] x, Quantity[] y) {
		if (x.length != y.length)
			throw new IllegalArgumentException("These array lengths don't match.");
		try {
			return interp(x, quadargmax(Math.max(0, left), Math.min(x.length, right), y));
		} catch (IndexOutOfBoundsException e) { // y is empty or all NaN
			return new Quantity(-1, y[0].getN());
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
		int i0 = Math.max(0, Math.min(x.length-2, (int) i));
		return (i0+1-i)*x[i0] + (i-i0)*x[i0+1];
	}
	
	/**
	 * take the floating-point index of an array using linear interpolation.
	 * @param x the array of values
	 * @param i the partial index
	 * @return x[i], more or less
	 */
	public static Quantity interp(double[] x, Quantity i) {
		if (i.value < 0 || i.value > x.length-1)
			throw new IndexOutOfBoundsException("Even partial indices have limits: "+i);
		int i0 = Math.max(0, Math.min(x.length-2, (int) i.value));
		return i.minus(i0).times(x[i0+1]).minus(i.minus(i0+1).times(x[i0]));
	}
	
	/**
	 * take the floating-point index of an array using linear interpolation.
	 * @param x the array of values
	 * @param i the partial index
	 * @return x[i], more or less
	 */
	public static Quantity interp(Quantity[] x, Quantity i) {
		if (i.value < 0 || i.value > x.length-1)
			throw new IndexOutOfBoundsException("Even partial indices have limits: "+i);
		int i0 = Math.max(0, Math.min(x.length-2, (int) i.value));
		return i.minus(i0).times(x[i0+1]).minus(i.minus(i0+1).times(x[i0]));
	}
	
	/**
	 * interpolate the value onto the given array.
	 * @param x0 the desired coordinate
	 * @param x the array of coordinates (must be unimodally increasing)
	 * @param y the array of values
	 * @return y(x0), more or less
	 */
	public static Quantity interp(Quantity x0, double[] x, Quantity[] y) {
		if (x0.value < x[0] || x0.value > x[x.length-1])
			throw new IndexOutOfBoundsException("Nope. Not doing extrapolation: "+x0);
		int l = 0, r = x.length;
		while (r - l > 1) {
			int m = (l + r)/2;
			if (x0.value < x[m])
				r = m;
			else
				l = m;
		}
		return y[l].times(x0.minus(x[r]).over(x[l] - x[r])).plus(y[r].times(x0.minus(x[l]).over(x[r] - x[l])));
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
		double[] dydx = new double[x.length-1];
		for (int i = 0; i < dydx.length; i ++)
			dydx[i] = (y[i+1] - y[i])/(x[i+1] - x[i]);
		return dydx;
	}
	
	/**
	 * find the second order finite difference derivative. for best results, x
	 * should be evenly spaced.
	 * @param x the x values
	 * @param y the corresponding y values
	 * @return the slope dy/dx at each point
	 */
	public static Quantity[] derivative(double[] x, Quantity[] y) {
		if (x.length != y.length)
			throw new IllegalArgumentException("Array lengths do not correspond.");
		if (x.length < 3)
			throw new IllegalArgumentException("I can't make inferences in these condicions!");
		Quantity[] dydx = new Quantity[x.length];
		for (int i = 0; i < y.length; i ++) {
			if (i == 0)
				dydx[i] = y[i].times(-3).plus(y[i+1].times(4)).plus(y[i+2].times(-1)).over(x[i+2] - x[i]);
			else if (i < y.length - 1)
				dydx[i] = y[i+1].minus(y[i-1]).over(x[i+1] - x[i-1]);
			else
				dydx[i] = y[i-2].times(-1).plus(y[i-1].times(4)).plus(y[i].times(-3)).over(x[i] - x[i-2]);
		}
		return dydx;
	}
	
	/**
	 * find the second order finite difference double derivative. for best results, x
	 * should be evenly spaced.
	 * @param x the x values
	 * @param y the corresponding y values
	 * @return the slope d^2y/dx^2 at each point
	 */
	public static Quantity[] secondDerivative(double[] x, Quantity[] y) {
		if (x.length != y.length)
			throw new IllegalArgumentException("Array lengths do not correspond.");
		if (x.length < 3)
			throw new IllegalArgumentException("I can't make inferences in these condicions!");
		Quantity[] dydx = new Quantity[x.length-1];
		for (int i = 0; i < y.length; i ++) {
			if (i == 0)
				dydx[i] = y[i].plus(y[i+1].times(-2)).plus(y[i+2]).over(Math.pow(x[i+1] - x[i-1], 2)/4.);
			else if (i < y.length - 1)
				dydx[i] = y[i+1].minus(y[i].times(2)).plus(y[i-1]).over(Math.pow(x[i+1] - x[i-1], 2)/4.);
			else
				dydx[i] = y[i-2].plus(y[i-1].times(-2)).plus(y[i]).over(Math.pow(x[i+1] - x[i-1], 2)/4.);
		}
		return dydx;
	}
	
	/**
	 * find the finite difference derivative. for best results, x should be evenly spaced.
	 * @param x the x values
	 * @param y the corresponding y values
	 * @param x the time at which to computer it
	 * @param Δx the time interval over which to compute it
	 * @return the slope dy/dx at each point
	 */
	public static Quantity derivative(double[] x, Quantity[] y, Quantity x0, double Δx) {
		if (x.length != y.length)
			throw new IllegalArgumentException("Array lengths do not correspond.");
		return interp(x0.plus(Δx/2.), x, y).minus(interp(x0 .minus(Δx/2.), x, y)).over(Δx);
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
	 * only works if the input spectrum has a higher resolution than the output spectrum :P
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
		
		double[][] zO = new double[yO.length-1][xO.length-1]; // resize the input array to match the output array
		for (int iI = 0; iI < yI.length-1; iI ++) {
			for (int jI = 0; jI < xI.length-1; jI ++) { // for each small pixel on the input spectrum
				double iO = (yI[iI] - yO[0])/(yO[1] - yO[0]); // find the big pixel of the scaled spectrum
				double jO = (xI[jI] - xO[0])/(xO[1] - xO[0]); // that contains the upper left corner
				int iOint = (int) Math.floor(iO);
				double iOmod = iO - iOint;
				int jOint = (int) Math.floor(jO);
				double jOmod = jO - jOint;
				double cU = Math.min(1, (1 - iOmod)*(yO[1] - yO[0])/(yI[iI+1] - yI[iI])); // find the fraction of it that is above the next pixel
				double cL = Math.min(1, (1 - jOmod)*(xO[1] - xO[0])/(xI[jI+1] - xI[jI])); // and left of the next pixel
				
				addIfInBounds(zO, iOint,   jOint,   zI[iI][jI]*cU*cL); // now add the contents of this spectrum
				addIfInBounds(zO, iOint,   jOint+1, zI[iI][jI]*cU*(1-cL)); // being careful to distribute them properly
				addIfInBounds(zO, iOint+1, jOint,   zI[iI][jI]*(1-cU)*cL); // (I used this convenience method because otherwise I would have to check all the bounds all the time)
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
	 * extract the values from an array of Quantities
	 * @return the value of each Quantity in the same order as before
	 */
	public static double[] modes(Quantity[] x) {
		double[] y = new double[x.length];
		for (int i = 0; i < x.length; i ++)
			y[i] = x[i].value;
		return y;
	}
	
	/**
	 * extract the errors from an array of Quantities
	 * @return the standard deviation of each Quantity in the same order as before
	 */
	public static double[] stds(Quantity[] x, double[][] covariance) {
		double[] y = new double[x.length];
		for (int i = 0; i < x.length; i ++)
			y[i] = Math.sqrt(x[i].variance(covariance));
		return y;
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
	 * copied from https://www.sanfoundry.com/java-program-find-inverse-matrix/
	 * @param a
	 * @return
	 */
	public static double[][] matinv(double[][] arr) {
		double[][] a = new double[arr.length][];
		for (int i = 0; i < arr.length; i ++) {
			if (arr[i].length != arr.length)
				throw new IllegalArgumentException("Only square matrices have inverses; not this "+arr.length+"×"+arr[i].length+" trash.");
			a[i] = arr[i].clone();
		}
		
		int n = a.length;
		double x[][] = new double[n][n];
		double b[][] = new double[n][n];
		int index[] = new int[n];
		for (int i = 0; i < n; ++i)
			b[i][i] = 1;

		// Transform the matrix into an upper triangle
		gaussian(a, index);

		// Update the matrix b[i][j] with the ratios stored
		for (int i = 0; i < n - 1; ++i)
			for (int j = i + 1; j < n; ++j)
				for (int k = 0; k < n; ++k)
					b[index[j]][k] -= a[index[j]][i] * b[index[i]][k];

		// Perform backward substitutions
		for (int i = 0; i < n; ++i) {
			x[n - 1][i] = b[index[n - 1]][i] / a[index[n - 1]][n - 1];
			for (int j = n - 2; j >= 0; --j) {
				x[j][i] = b[index[j]][i];
				for (int k = j + 1; k < n; ++k) {
					x[j][i] -= a[index[j]][k] * x[k][i];
				}
				x[j][i] /= a[index[j]][j];
			}
		}
		return x;
	}
	
	/**
	 * Method to carry out the partial-pivoting Gaussian
	 * elimination. Here index[] stores pivoting order.
	 * @param a
	 * @param index
	 */
	private static void gaussian(double a[][], int index[]) {
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

			// Interchange rows according to the pivoting order
			int itmp = index[j];
			index[j] = index[k];
			index[k] = itmp;
			for (int i = j + 1; i < n; ++i) {
				double pj = a[index[i]][j] / a[index[j]][j];

				// Record pivoting ratios below the diagonal
				a[index[i]][j] = pj;

				// Modify other elements accordingly
				for (int l = j + 1; l < n; ++l)
					a[index[i]][l] -= pj * a[index[j]][l];
			}
		}
	}
	
	
	/**
	 * a poor person's pseudoinverse. It's like a regular inverse, but if a particular diagonal
	 * value is zero, then it removes that dimension before inverting, and then puts NaNs back
	 * in where they were.
	 * @param arr
	 * @return
	 */
	public static double[][] pseudoinv(double[][] arr) {
		int n = 0;
		boolean[] useful = new boolean[arr.length];
		for (int i = 0; i < arr.length; i ++) {
			useful[i] = (Double.isFinite(arr[i][i]) && arr[i][i] != 0);
			if (useful[i])  n ++;
		}
		double[][] a = new double[n][n];
		int k = 0;
		for (int i = 0; i < arr.length; i ++) {
			if (useful[i]) {
				int l = 0;
				for (int j = 0; j < arr[i].length; j ++) {
					if (useful[j]) {
						a[k][l] = arr[i][j];
						l ++;
					}
				}
				k ++;
			}
		}
		double[][] b = matinv(a);
		double[][] c = new double[arr.length][arr.length];
		k = 0;
		for (int i = 0; i < arr.length; i ++) {
			if (useful[i]) {
				int l = 0;
				for (int j = 0; j < arr[i].length; j ++) {
					if (useful[j]) {
						c[i][j] = b[k][l];
						l ++;
					}
					else {
						c[i][j] = 0;
					}
				}
				k ++;
			}
			else {
				for (int j = 0; j < arr[i].length; j ++)
					c[i][j] = 0;
			}
		}
		return c;
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
	 * Legendre polynomial of degree n
	 * @param n the order of the polynomial
	 * @param z the cosine of the angle at which this is evaluated
	 * @return P_l(z)
	 */
	public static double legendre(int n, double z) {
		if (n == 0)
			return 1;
		else if (n == 1)
			return z;
		else if (n == 2)
			return (3*z*z - 1)/2.;
		else if (n == 3)
			return (5*z*z - 3)*z/2.;
		else if (n == 4)
			return ((35*z*z - 30)*z*z + 3)/8.;
		else if (n == 5)
			return ((63*z*z - 70)*z*z + 15)*z/8.;
		else if (n == 6)
			return (((231*z*z - 315)*z*z + 105)*z*z - 5)/16.;
		else if (n == 7)
			return (((429*z*z - 693)*z*z + 315)*z*z - 35)*z/16.;
		else if (n == 8)
			return ((((6435*z*z - 12012)*z*z + 6930)*z*z - 1260)*z*z + 35)/128.;
		else if (n == 9)
			return ((((12155*z*z - 25740)*z*z + 18018)*z*z - 4620)*z*z + 315)*z/128.;
		else
			throw new IllegalArgumentException("I don't know Legendre polynomials that high.");
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
		 * instantiate a new function given x and y data in columns, and assuming x values are
		 * all equally spaced
		 * @param data array of {x, y}
		 * @param resolution the number of x intervals
		 */
		public DiscreteFunction(double[][] data, int resolution) {
			if (resolution != data.length-1)
				throw new IllegalArgumentException("this resolution is a lie");
			for (int i = 1; i < data.length; i ++)
				if (data[i][0] < data[i-1][0])
					throw new IllegalArgumentException("x must be monotonically increasing.");
			
			this.X = new double[data.length];
			this.Y = new double[data.length];
			for (int i = 0; i < data.length; i ++) {
				X[i] = data[i][0];
				Y[i] = data[i][1];
			}
			this.resolution = resolution;
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
			if (resolution != x.length-1)
				throw new IllegalArgumentException("this resolution is a lie");
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
			else if (x >= X[X.length-1]) // or highest values, depending on which is appropriate
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
		 * it's a function. evaluate it. if this function's x values are equally spaced, this
		 * can be run in O(1) time. otherwise, it will take O(log(n)).
		 * @param x
		 * @return f(x)
		 */
		public Quantity evaluate(Quantity x) {
			int i; // we will linearly interpolate x from (X[i], X[i+1]) onto (Y[i], Y[i+1]).
			if (x.value < X[0]) // if it's out of bounds, we will extrapolate from the lowest values
				i = 0;
			else if (x.value > X[X.length-1]) // or highest values, depending on which is appropriate
				i = X.length-2;
			else if (this.resolution > 0) // nonzero resolution means we can find i itself with linear interpolation
				i = (int)((x.value - X[0])/(X[resolution] - X[0])*resolution); // linearly interpolate x from X to i
			else { // otherwise, we'll need a binary search
				int min = 0, max = X.length; // you know about binary searches, right?
				i = (min + max)/2;
				while (max - min > 1) { // I probably don't need to explain this.
					if (X[i] < x.value)
						min = i;
					else
						max = i;
					i = (min + max)/2;
				}
			}
			return x.minus(X[i]).times((Y[i+1] - Y[i])/(X[i+1] - X[i])).plus(Y[i]); // linearly interpolate x from X[i] to Y[i]
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
	
	
	/**
	 * A value that tracks its gradient in parameter space for the purpose of error bar
	 * determination.
	 * @author Justin Kunimune
	 */
	public static class Quantity {
		public final double value;
		public final Vector gradient;
		
		public Quantity(double value, int n) {
			this(value, new double[n]);
		}
		
		public Quantity(double value, double[] gradient) {
			this(value, new Vector(gradient));
		}

		public Quantity(double value, Vector gradient) {
			if (Double.isNaN(value))
				throw new IllegalArgumentException("I did not account for this!");
			this.value = value;
			this.gradient = gradient;
		}
		
		public double variance(double[][] covariance) {
			double variance = this.gradient.dot(new Matrix(covariance).times(this.gradient));
			if (variance < 0) { // if it doesn't work
				double[][] newCovariance = new double[covariance.length][covariance.length];
				for (int i = 0; i < covariance.length; i ++) {
					if (covariance[i][i] < 0)
						return variance; // first check that the diagonal terms are positive (they should be)
					for (int j = 0; j < covariance.length; j ++) {
						if (i == j)  newCovariance[i][j] = covariance[i][j]; // then halving the off-diagonal terms and try again
						else         newCovariance[i][j] = covariance[i][j]/2;
					}
				}
				return this.variance(newCovariance);
			}
			else {
				return variance;
			}
		}
		
		public Quantity plus(double constant) {
			return new Quantity(this.value + constant, this.gradient);
		}
		
		public Quantity plus(Quantity that) {
			return new Quantity(this.value + that.value, this.gradient.plus(that.gradient));
		}
		
		public Quantity minus(double constant) {
			return new Quantity(this.value - constant, this.gradient);
		}
		
		public Quantity minus(Quantity that) {
			return this.plus(that.times(-1));
		}
		
		public Quantity times(double constant) {
			return new Quantity(this.value*constant, this.gradient.times(constant));
		}
		
		public Quantity times(Quantity that) {
			return new Quantity(this.value*that.value,
					this.gradient.times(that.value).plus(that.gradient.times(this.value)));
		}
		
		public Quantity over(double constant) {
			return this.times(1/constant);
		}
		
		public Quantity over(Quantity that) {
			return new Quantity(this.value/that.value,
					this.gradient.times(that.value).minus(that.gradient.times(this.value)).times(
							Math.pow(that.value, -2)));
		}
		
		public Quantity pow(double exponent) {
			return new Quantity(Math.pow(this.value, exponent),
					this.gradient.times(exponent*Math.pow(this.value, exponent - 1)));
		}
		
		public Quantity sqrt() {
			return this.pow(1/2.);
		}
		
		public Quantity abs() {
			if (this.value < 0)
				return this.times(-1);
			else
				return this;
		}
		
		public Quantity mod(double divisor) {
			return new Quantity(this.value%divisor, this.gradient);
		}
		
		/**
		 * @return the number of variable dimensions in which this exists
		 */
		public int getN() {
			return this.gradient.getN();
		}
		
		public String toString(double[][] covariance) {
			return String.format("%8.6g \u00B1 %8.3g", this.value, Math.sqrt(this.variance(covariance)));
		}
	}
	
	
	private static class Matrix {
		private final double[][] values;
		
		public Matrix(double[][] values) {
			this.values = values;
		}
		
		public Vector times(Vector v) {
			if (v.getN() != this.getM())
				throw new IllegalArgumentException("the dimensions don't match.");
			double[] product = new double[this.getN()];
			for (int i = 0; i < this.getN(); i ++)
				for (int j = 0; j < this.getM(); j ++)
					if (this.values[i][j] != 0 && v.values[j] != 0) // 0s override Infs and NaNs in this product
						product[i] += this.values[i][j]*v.values[j];
			return new Vector(product);
		}
		
		public int getN() {
			return this.values.length;
		}
		
		public int getM() {
			return this.values[0].length;
		}
	}
	
	private static class Vector {
		private final double[] values;
		
		public Vector(double[] values) {
			this.values = values;
		}
		
		public Vector plus(Vector that) {
			double[] sum = new double[this.getN()];
			for (int i = 0; i < this.getN(); i ++)
				sum[i] = this.values[i] + that.values[i];
			return new Vector(sum);
		}
		
		public Vector minus(Vector that) {
			return this.plus(that.times(-1));
		}
		
		public Vector times(double scalar) {
			double[] product = new double[this.getN()];
			for (int i = 0; i < this.getN(); i ++)
				product[i] = this.values[i]*scalar;
			return new Vector(product);
		}
		
		public double dot(Vector that) {
			if (this.getN() != that.getN())
				throw new IllegalArgumentException("the dimensions don't match.");
			double product = 0;
			for (int i = 0; i < this.getN(); i ++)
				if (this.values[i] != 0 && that.values[i] != 0)
					product += this.values[i]*that.values[i];
			return product;
		}

		public int getN() {
			return values.length;
		}
	}
	
	
	public static final void main(String[] args) {
		double[][] cov = {{1, 0}, {0, 1}};
		Quantity x = new Quantity(5, new double[] {1, 0});
		Quantity y = new Quantity(12, new double[] {0, 1});
		System.out.println(x.toString(cov));
		System.out.println(y.toString(cov));
		System.out.println(x.plus(y).toString(cov));
		System.out.println(x.minus(y).toString(cov));
		System.out.println(x.times(y).toString(cov));
		System.out.println(x.over(y).toString(cov));
		System.out.println(x.mod(4).toString(cov));
	}
	
}
