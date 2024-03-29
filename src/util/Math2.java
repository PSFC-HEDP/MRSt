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
package util;

import java.util.Arrays;
import java.util.Locale;
import java.util.Random;
import java.util.function.Function;

/**
 * a file with some useful numerical analysis stuff.
 * 
 * @author Justin Kunimune
 */
public class Math2 {

	/**
	 * draw a boolean from a Bernoulli distribution.
	 * @param p the probability of true
	 * @return the number
	 */
	public static boolean bernoulli(double p) {
		return bernoulli(p, Math.random());
	}

	/**
	 * draw a boolean from a Bernoulli distribution.
	 * @param p the probability of true
	 * @param random the rng to use
	 * @return the number
	 */
	public static boolean bernoulli(double p, Random random) {
		return bernoulli(p, random.nextDouble());
	}

	/**
	 * draw a boolean from a Bernoulli distribution using the given random number.
	 * @param p the probability of true
	 * @param u a number randomly distributed in [0, 1)
	 * @return the number
	 */
	private static boolean bernoulli(double p, double u) {
		return u < p;
	}

	/**
	 * draw a number from a Gaussian distribution.
	 * @param μ mean
	 * @param σ standard deviation
	 * @param random the rng to use
	 * @return the number
	 */
	public static double normal(double μ, double σ, Random random) {
		return μ + σ*random.nextGaussian();
	}

	/**
	 * draw a number from a Poisson distribution.
	 * @param λ expected number
	 * @param random the rng to use
	 * @return the number
	 */
	public static int poisson(double λ, Random random) {
		if (λ < 20)
			return poisson(λ, random.nextDouble());
		else if (λ < 2e+9)
			return (int) Math.max(0., Math.round(normal(λ, Math.sqrt(λ), random)));
		else
			throw new IllegalArgumentException("this method will hit runoff issues if you go so hi ("+λ+").");
	}

	/**
	 * draw a number from a Poisson distribution using the given random number.
	 * @param λ expected number
	 * @param u a number randomly distributed in [0, 1)
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

	/**
	 * draw a number from an exponential distribution.
	 * @param μ mean
	 * @return the number
	 */
	public static double exponential(double μ) {
		return exponential(μ, Math.random());
	}

	/**
	 * draw a number from an exponential distribution.
	 * @param μ mean
	 * @param random the rng to use
	 * @return the number
	 */
	public static double exponential(double μ, Random random) {
		return exponential(μ, random.nextDouble());
	}

	/**
	 * draw a number from an exponential distribution using the given random number.
	 * @param λ mean
	 * @param u a number randomly distributed in [0, 1)
	 * @return the number
	 */
	private static double exponential(double λ, double u) {
		return -λ*Math.log(1 - u);
	}

	/**
	 * draw a number from an erlang distribution.
	 * @param k number of exponential distributions to sum
	 * @param μ mean of each individual exponential distribution
	 * @param random the rng to use
	 * @return the number
	 */
	public static double erlang(int k, double μ, Random random) {
		if (k < 20) {
			double u = 1;
			for (int i = 0; i < k; i ++)
				u *= 1 - random.nextDouble();
			return exponential(μ, 1 - u);
		}
		else {
			return Math.max(0., normal(k*μ, Math.sqrt(k)*μ, random));
		}
	}

	public static double gamma(double a, double b, Random random) {
		if (a <= 0) {
			throw new IllegalArgumentException("a must > 0");
		}
		if (a < 20) {
			double μ = a/b, σ = Math.sqrt(a)/b;
			double[] x = new double[40];
			double[] y = new double[x.length-1];
			for (int i = 0; i < x.length; i ++)
				x[i] = Math.max(0., μ - 4.*σ) + 8.*σ/(x.length-1)*i;
			for (int i = 0; i < y.length; i ++) {
				double xM = (x[i] + x[i+1])/2.;
				y[i] = Math.exp((a-1)*Math.log(xM) - b*xM - ((a-1)*Math.log(μ) - b*μ));
			}
			return drawFromProbabilityDistribution(x, y, random);
		}
		else {
			return Math.max(0., normal(a/b, Math.sqrt(a)/b, random));
		}
	}

	public static int binomial(int n, double p, Random random) {
		if (p < 0 || p > 1)
			throw new IllegalArgumentException("p must be in [0, 1] but you passd "+p);
		else if (p == 0)
			return 0;
		else if (p == 1)
			return n;
		if (n < 0)
			throw new IllegalArgumentException("n must > 0 but you passd "+n);
		else if (n == 0)
			return 0;
		else if (p*n < 20) {
			double[] P = new double[n + 1];
			double logp = Math.log(p);
			double logq = Math.log(1 - p);
			for (int i = 0; i <= Math.min(n, 40); i ++)
				P[i] = Math.exp(Math2.logChoose(n, i) + i*logp + (n - i)*logq);
			double total = Math2.sum(P);
			double u = random.nextDouble()*total;
			for (int i = 0; i <= n; i ++) {
				if (u < P[i])
					return i;
				else
					u -= P[i];
			}
			throw new IllegalArgumentException("math is broken: "+n+","+p+":"+ Arrays.toString(P)+", "+total);
		}
		else if ((1 - p)*n < 20) {
			return n - binomial(n, 1 - p, random);
		}
		else {
			return (int) Math.round(Math.max(0, Math.min(n,
					normal(n*p, Math.sqrt(n*p*(1 - p)), random))));
		}
	}

	public static double drawFromProbabilityDistribution(
		  double[] x, double[] pdf, Random random) {
		double sum = 0;
		for (int i = 0; i < pdf.length; i ++)
			sum += pdf[i]*(x[i+1] - x[i]);
		double u = random.nextDouble();
		for (int i = 0; i < pdf.length; i ++) {
			double thresh = pdf[i]*(x[i + 1] - x[i])/sum;
			if (u <= thresh)
				return x[i] + (x[i+1] - x[i])/thresh*u;
			else
				u -= thresh;
		}
		throw new IllegalArgumentException("math is broken: "+ Arrays.toString(x)+", "+Arrays.toString(pdf));
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
		return average(y, f, 0, y.length);
	}
	
	/**
	 * compute the averaged value in the given interval
	 * @param y the values
	 * @param f the weights
	 * @param left the starting index (inclusive)
	 * @param rite the ending index (exclusive)
	 */
	public static Quantity average(Quantity[] y, Quantity[] f, int left, int rite) {
		if (y.length != f.length)
			throw new IllegalArgumentException("The array lengths don't match.");
		Quantity p0 = new Quantity(0, y[0].getN());
		Quantity p1 = new Quantity(0, y[0].getN());
		for (int i = left; i < rite; i ++) {
			p0 = p0.plus(f[i]);
			p1 = p1.plus(f[i].times(y[i]));
		}
		return p1.over(p0);
	}

	/**
	 * find the full-width at half-maximum of a distribucion. if it is very noisy, this will
	 * underestimate the width.
	 * @param x the bin edges
	 * @param y the number in each bin
	 * @return the full-width at half-maximum
	 */
	public static double fwhm(double[] x, double[] y) {
		Quantity[] y_q = new Quantity[y.length];
		for (int i = 0; i < y.length; i ++)
			y_q[i] = new Quantity(y[i], 0);
		Quantity width = fwhm(x, y_q);
		if (width == null)
			return Double.POSITIVE_INFINITY;
		else
			return width.value;
	}

	/**
	 * find the full-width at half-maximum of a distribucion. if it is very noisy, this will
	 * underestimate the width.
	 * @param x the bin edges
	 * @param y the number in each bin
	 * @return the full-width at half-maximum
	 */
	public static Quantity fwhm(double[] x, Quantity[] y) {
		if (x.length != y.length + 1)
			throw new IllegalArgumentException("the inputs must have matching lengths.");
		x = Math2.binCenters(x);
		int max = argmax(y);
		Quantity xR = null;
		for (int i = max + 1; i < x.length; i ++) {
			if (y[i].value < y[max].value/2.) {
				xR = interp(y[max].over(2.), y[i-1], y[i], x[i-1], x[i]);
				break;
			}
		}
		Quantity xL = null;
		for (int i = max; i >= 1; i --) {
			if (y[i-1].value < y[max].value/2.) {
				xL = interp(y[max].over(2.), y[i-1], y[i], x[i-1], x[i]);
				break;
			}
		}
		if (xR == null || xL == null)
			return new Quantity(Double.POSITIVE_INFINITY, y[0].getN());
		else {
			if (xR.minus(xL).value < (x[1] - x[0])*4)
				System.err.println("Warning: this FWHM estimate isn't binned enuff.");
			return xR.minus(xL);
		}
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
	 * @param x the array of points
	 * @return the standard deviation
	 */
	public static double std(double[] x) {
		double mean = 0;
		double meanSqr = 0;
		for (double v: x) {
			mean += v/x.length;
			meanSqr += Math.pow(v, 2)/x.length;
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

	public static double[] sum(double[][] arr, int axis) {
		if (axis == 0) {
			double[] s = new double[arr[0].length];
			for (double[] row: arr)
				for (int j = 0; j < row.length; j++)
					s[j] += row[j];
			return s;
		}
		else if (axis == 1) {
			double[] output = new double[arr.length];
			for (int i = 0; i < arr.length; i ++)
				output[i] = Math2.sum(arr[i]);
			return output;
		}
		else {
			throw new IllegalArgumentException(Integer.toString(axis));
		}
	}

	public static double mean(double[] arr) {
		return sum(arr)/arr.length;
	}

	/**
	 * compute the sum of a normalized histogram
	 * @param x the bin edges corresponding to the outer index
	 * @param y the density in each bin
	 * @return the total value
	 */
	public static double integral(double[] x, double[] y) {
		if (x.length != y.length + 1)
			throw new IllegalArgumentException("the array lengths do not correspond");
		double sum = 0;
		for (int i = 0; i < x.length - 1; i ++)
				sum += y[i]*(x[i + 1] - x[i]);
		return sum;
	}

	/**
	 * compute the sum of a normalized histogram
	 * @param x the bin edges corresponding to the outer index
	 * @param y the bin edges corresponding to the inner index
	 * @param z the density in each bin
	 * @return the total value
	 */
	public static double iintegral(double[] x, double[] y, double[][] z) {
		if (x.length != z.length + 1 || y.length != z[0].length + 1)
			throw new IllegalArgumentException("for a "+z.length+"x"+z[0].length+" array, the bin edges cannot be "+x.length+" and "+y.length);
		double sum = 0;
		for (int i = 0; i < x.length - 1; i ++)
			for (int j = 0; j < y.length - 1; j ++)
				sum += z[i][j]*(x[i + 1] - x[i])*(y[j + 1] - y[j]);
		return sum;
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
	
	public static double sqr(double[] v) {
		double s = 0;
		for (double x: v)
			s += Math.pow(x, 2);
		return s;
	}
	
	public static int lastIndexBefore(double level, double[] v, int start) {
		int l = start;
		while (l-1 >= 0 && v[l-1] > level)
			l --;
		return l;
	}
	
	public static int firstIndexAfter(double level, double[] v, int start) {
		int r = start;
		while (r < v.length && v[r] > level)
			r ++;
		return r;
	}

	public static boolean[][] nonzero(double[][] values) {
		boolean[][] nonzero = new boolean[values.length][values[0].length];
		for (int i = 0; i < values.length; i ++)
			for (int j = 0; j < values[i].length; j ++)
				nonzero[i][j] = values[i][j] != 0;
		return nonzero;
	}

	public static double[][] deepCopy(double[][] orig) {
		double[][] copy = new double[orig.length][];
		for (int i = 0; i < copy.length; i ++)
			copy[i] = Arrays.copyOf(orig[i], orig[i].length);
		return copy;
	}

	public static double min(double[] arr) {
		double min = Double.POSITIVE_INFINITY;
		for (double x: arr)
			if (x < min)
				min = x;
		return min;
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
	 * find the last index of the highest value
	 * @param x the array of values
	 * @return i such that x[i] >= x[j] for all j
	 */
	public static int argmax(Quantity[] x) {
		int argmax = -1;
		for (int i = 0; i < x.length; i ++)
			if (!Double.isNaN(x[i].value) && (argmax == -1 || x[i].value > x[argmax].value))
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
			throw new IndexOutOfBoundsException("Even partial indices have limits: "+i+"/"+x.length);
		int i0 = Math.max(0, Math.min(x.length-2, (int) i));
		return (i0+1-i)*x[i0] + (i-i0)*x[i0+1];
	}

	/**
	 * interpolate a value onto a line
	 */
	public static double interp(double x, double x1, double x2, double y1, double y2) {
		return y1 + (x - x1)/(x2 - x1)*(y2 - y1);
	}

	/**
	 * interpolate a value onto a line
	 */
	public static Quantity interp(Quantity x, Quantity x1, Quantity x2, double y1, double y2) {
		return x.minus(x1).over(x2.minus(x1)).times(y2 - y1).plus(y1);
	}

	/**
	 * take the floating-point index of an array using linear interpolation.
	 * @param x the array of values
	 * @param i the partial index
	 * @return x[i], more or less
	 */
	public static Quantity interp(double[] x, Quantity i) {
		Quantity[] xQ = new Quantity[x.length];
		for (int j = 0; j < x.length; j ++)
			xQ[j] = new Quantity(x[j], i.getN());
		return interp(xQ, i);
	}
	
	/**
	 * take the floating-point index of an array using linear interpolation.
	 * @param x the array of values
	 * @param i the partial index
	 * @return x[i], more or less
	 */
	public static Quantity interp(Quantity[] x, Quantity i) {
		if (i.value < 0 || i.value > x.length-1)
			throw new IndexOutOfBoundsException("Even partial indices have limits: "+i.value+"/"+x.length);
		if (x.length == 1)
			return x[0];
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
	public static double interp(double x0, double[] x, double[] y) {
		Quantity[] yQ = new Quantity[y.length];
		for (int j = 0; j < y.length; j ++)
			yQ[j] = new Quantity(y[j], 0);
		return interp(new Quantity(x0, 0), x, yQ).value;
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
	 * take the floating-point index of an array using cubic interpolation.
	 * @param x the array of values
	 * @param i the partial index
	 * @return x[i], more or less
	 */
	public static Quantity quadInterp(Quantity[] x, Quantity i) {
		if (i.value < 0 || i.value > x.length-1)
			throw new IndexOutOfBoundsException("Even partial indices have limits: "+i.value+"/"+x.length);
		if (x.length == 1)
			return x[0];
		else if (x.length < 4)
			throw new UnsupportedOperationException("I haven't implemented this and don't want to.");
		int i0 = Math.max(1, Math.min(x.length-3, (int) i.value));
		Quantity xA = x[i0-1], xB = x[i0], xC = x[i0+1], xD = x[i0+2];
		Quantity δA = i.minus(i0 - 1), δB = i.minus(i0), δC = i.minus(i0 + 1).times(-1), δD = i.minus(i0 + 2).times(-1);
		return xB.times(δA).times(δC).times(δD).times(3).plus(xC.times(δA).times(δB).times(δD).times(3)).minus(xA.times(δB).times(δC).times(δD)).minus(xD.times(δA).times(δB).times(δC)).over(6);
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
		Quantity[] d2ydx2 = new Quantity[x.length];
		for (int i = 1; i < x.length - 1; i ++)
			d2ydx2[i] = y[i+1].minus(y[i].times(2)).plus(y[i-1]).over(Math.pow(x[i+1] - x[i-1], 2)/4.);
		d2ydx2[0] = d2ydx2[1];
		d2ydx2[x.length-1] = d2ydx2[x.length-2];
		return d2ydx2;
	}
	
	private static Quantity[][] getLocalDataProperties(double[] x, Quantity[] y, Quantity x0, double Δx) {
		if (x.length != y.length)
			throw new IllegalArgumentException("Array lengths do not correspond.");
		int n = y[0].getN();
		Quantity xL = x0.minus(Δx/2), xR = x0.plus(Δx/2);
		if (xL.value < x[0])
			xL = new Quantity(x[0], n);
		if (xR.value > x[x.length-1])
			xR = new Quantity(x[x.length-1], n);
		
		int s = 0;
		for (int i = 0; i < x.length; i ++)
			if (x[i] > xL.value && x[i] < xR.value)
				s ++;
		s += 2;
		
		Quantity[] xData = new Quantity[s];
		int l = 1;
		for (int i = 0; i < x.length; i ++) {
			if (x[i] > xL.value && x[i] < xR.value) {
				xData[l] = new Quantity(x[i], n);
				l ++;
			}
		}
		xData[0] = xL;
		xData[s-1] = xR;
		
		Quantity[] wate = new Quantity[s];
		wate[0] = xData[1].minus(xData[0]);
		for (int i = 1; i < s-1; i ++)
			wate[i] = xData[i+1].minus(xData[i-1]);
		wate[s-1] = xData[s-1].minus(xData[s-2]);
		
		Quantity[] yData = new Quantity[s];
		for (int i = 0; i < s; i ++)
			yData[i] = interp(xData[i], x, y);
		Quantity[][] sums = new Quantity[5][2];
		for (int j = 0 ; j < sums.length; j ++)
			for (int k = 0; k < sums[j].length; k ++)
				sums[j][k] = new Quantity(0, n);
		for (int i = 0; i < s; i ++)
			for (int j = 0; j < sums.length; j ++)
				for (int k = 0; k < sums[j].length; k ++)
					sums[j][k] = sums[j][k].plus(xData[i].pow(j).times(yData[i].pow(k)).times(wate[i]));
		System.out.println("looking between "+(x0.value - Δx/2)+" and "+(x0.value + Δx/2));
		System.out.println(s);
		
		Quantity[][] moments = new Quantity[5][2];
		for (int j = 0; j < moments.length; j ++)
			for (int k = 0; k < moments[j].length; k ++)
				moments[j][k] = sums[j][k].over(sums[0][0]);
		return moments;
	}
	
	/**
	 * fit to a parabola and find the nth derivative.  x must be evenly spaced.
	 */
	public static Quantity derivative(double[] x, Quantity[] y, Quantity i0, double Δx, int n) {
		if (x.length < 2)
			return new Quantity(Double.NaN, y[0].getN());
		double dx = x[1] - x[0];
		Quantity x0 = Math2.interp(x, i0);
		Quantity[] weights = new Quantity[x.length];
		for (int i = 0; i < x.length; i ++) {
			if (x[i] <= x0.minus(Δx/2 + dx/2).value)
				weights[i] = new Quantity(0, x0.getN());
			else if (x[i] <= x0.minus(Δx/2 - dx/2).value)
				weights[i] = x0.minus(Δx/2 + dx/2).minus(x[i]).over(-dx);
			else if (x[i] <= x0.plus(Δx/2 - dx/2).value)
				weights[i] = new Quantity(1, x0.getN());
			else if (x[i] <= x0.plus(Δx/2 + dx/2).value)
				weights[i] = x0.plus(Δx/2 + dx/2).minus(x[i]).over(dx);
			else
				weights[i] = new Quantity(0, x0.getN());
		}
		
		double[] xMoments = new double[5];
		Quantity[] yMoments = new Quantity[3];
		for (int j = 0; j < 3; j ++)
				yMoments[j] = new Quantity(0, x0.getN());
		for (int i = 0; i < x.length; i ++) {
			if (weights[i].value > 0) {
				for (int j = 0; j < 5; j ++)
					xMoments[j] = xMoments[j] +
							weights[i].value * Math.pow(x[i], j);
				for (int j = 0; j < 3; j ++)
					yMoments[j] = yMoments[j].plus(
							weights[i].times(y[i]).times(Math.pow(x[i], j)));
			}
		}
		
		if (n == 1) {
			return yMoments[0].times(xMoments[1]).minus(yMoments[1].times(xMoments[0])).over(
					xMoments[1]*xMoments[1] - xMoments[2]*xMoments[0]);
		}
		if (n == 2) {
			double[][] mat = new double[3][3];
			for (int i = 0; i < 3; i ++)
				for (int j = 0; j < 3; j ++)
					mat[i][j] = xMoments[2 + i - j];
			double[][] matInv = matinv(mat);
			return yMoments[0].times(matInv[0][0]).plus(
					yMoments[1].times(matInv[0][1])).plus(
					yMoments[2].times(matInv[0][2])).times(2);
		}
		else
			throw new IllegalArgumentException("I don't do that derivative.");
	}
	
	/**
	 * return the index of the pair of bin edges in an evenly spaced array that contains
	 * the value
	 * @return int in the range [0, bins.length-1), or -1 if it's out of range
	 */
	public static int bin(double value, double[] binEdges) {
		if (Double.isNaN(value))
			return -1;
		int bin = (int)((value - binEdges[0])/(binEdges[binEdges.length-1] - binEdges[0])*(binEdges.length-1));
		return (bin >= 0 && bin < binEdges.length-1) ? bin : -1;
	}
	

	public static double[] collum(double[][] matrix, int collumIndex) {
		double[] collum = new double[matrix.length];
		for (int i = 0; i < matrix.length; i ++) {
			if (collumIndex >= matrix[i].length)
				throw new IllegalArgumentException("the given matrix does not have enuff collums");
			collum[i] = matrix[i][collumIndex];
		}
		return collum;
	}


	public static double[] full(double value, int size) {
		double[] array = new double[size];
		for (int i = 0; i < size; i ++)
			array[i] = value;
		return array;
	}


	public static Quantity[] full(Quantity value, int size) {
		Quantity[] array = new Quantity[size];
		for (int i = 0; i < size; i ++)
			array[i] = value;
		return array;
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

	public static double[][] transpose(double[] collum) {
		double[][] output = new double[collum.length][];
		for (int i = 0; i < collum.length; i ++)
			output[i] = new double[] {collum[i]};
		return output;
	}

	public static double[] binCenters(double[] edges) {
		double[] centers = new double[edges.length - 1];
		for (int i = 0; i < edges.length - 1; i ++)
			centers[i] = (edges[i] + edges[i + 1])/2;
		return centers;
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
	 * do a Runge-Kutta 4-5 integral to get the final value of y after some interval
	 * @param f dy/dt as a function of y
	 * @param Δt the amount of time to let y ject
	 * @param y0 the inicial value of y
	 * @param numSteps the number of steps to use
	 * @return the final value of y
	 */
	public static double odeSolve(DiscreteFunction f, double Δt, double y0, int numSteps) {
		final double dt = Δt/numSteps;
		double y = y0;
		for (int i = 0; i < numSteps; i ++) {
			double k1 = f.evaluate(y);
			double k2 = f.evaluate(y + k1/2.*dt);
			double k3 = f.evaluate(y + k2/2.*dt);
			double k4 = f.evaluate(y + k3*dt);
			y = y + (k1 + 2*k2 + 2*k3 + k4)/6.*dt;
		}
		return y;
	}
	
	/**
	 * a simple convenience method to avoid excessive if statements
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
	 * do the convolution of a vector with a 1d kernel.  b should be odd in length
	 * because I don't know how to center it otherwise, and symmetric because I
	 * don't know which way it's going to be oriented here.
	 * @param a vector
	 * @param b vector
	 * @return vector with the same length as a
	 */
	public static double[] convolve(double[] a, double[] b) {
		assert b.length%2 == 1;
		double[] c = new double[a.length];
		for (int i = 0; i < a.length; i ++)
			for (int j = 0; j < b.length; j ++)
				if (i + j - b.length/2 >= 0 && i + j - b.length/2 < c.length)
					c[i + j - b.length/2] += a[i]*b[j];
		return c;
	}

	/**
	 * multiply a vector by a matrix
	 * @param A matrix
	 * @param v vector
	 * @return A.v vector
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
	 * do a quick pass thru all of the 2x2 submatrices of this symmetric matrix
	 * to make sure they have nonnegative determinants, and alter the nondiagonal
	 * elements if they don't.
	 */
	public static void coercePositiveSemidefinite(double[][] A) {
		for (double[] row: A)
			if (row.length != A.length)
				throw new IllegalArgumentException("this method only works with square matrices.");

		for (int i = 0; i < A.length; i ++)
			if (A[i][i] < 0)
				A[i][i] = 0;

		for (int i = 0; i < A.length; i ++)
			for (int j = i+1; j < A.length; j ++)
				if (A[i][j]*A[j][i] > A[i][i]*A[j][j])
					A[i][j] = A[j][i] = Math.signum(A[i][j])*Math.sqrt(A[i][i]*A[j][j]); // enforce positive semidefiniteness
	}

	/**
	 * average this matrix with it transpose.
	 */
	public static void coerceSymmetric(double[][] A) {
		for (double[] row: A)
			if (row.length != A.length)
				throw new IllegalArgumentException("this method only works with square matrices.");

		for (int i = 0; i < A.length; i ++)
			for (int j = i+1; j < A.length; j ++)
				A[i][j] = A[j][i] = (A[i][j] + A[j][i])/2;
	}

	/**
	 * copied from <a href="https://www.sanfoundry.com/java-program-find-inverse-matrix/">sanfoundry.com</a>
	 */
	public static double[][] matinv(double[][] arr) {
		double[][] a = new double[arr.length][];
		for (int i = 0; i < arr.length; i ++) {
			if (arr[i].length != arr.length)
				throw new IllegalArgumentException("Only square matrices have inverses; not this "+arr.length+"×"+arr[i].length+" trash.");
			a[i] = arr[i].clone();
		}
		
		int n = a.length;
		double[][] x = new double[n][n];
		double[][] b = new double[n][n];
		int[] index = new int[n];
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
	 */
	private static void gaussian(double[][] a, int[] index) {
		int n = index.length;
		double[] c = new double[n];

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
	 * the log of the binomial coefficient function
	 * @param n the number of possibilities
	 * @param k the number in the selected set
	 * @return the number of combinations
	 */
	public static double logChoose(int n, int k) {
		double C = 1;
		for (int i = 1; i <= k; i ++)
			C += Math.log((double)(n - k + i)/i);
		return C;
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
	 * calculate the hessian matrix of a function using finite differences
	 * @param function the function to differentiate
	 * @param x0 the point at which to differentiate
	 * @param dx the step size of each dimension
	 * @return the symmetrick hessian matrix
	 */
	public static double[][] hessian(Function<double[], Double> function,
									 double[] x0, double[] dx) {
		if (x0.length != dx.length)
			throw new IllegalArgumentException("these arrays are supposed to have the same length.");

		double c = function.apply(x0); // start by getting the actual value
		assert Double.isFinite(c);

		double[] step = new double[x0.length]; // and the values in all basis directions
		for (int i = 0; i < step.length; i ++) {
			double[] xR = Arrays.copyOf(x0, x0.length);
			xR[i] += dx[i];
			step[i] = function.apply(xR);
			assert Double.isFinite(step[i]);

//			double min = Double.POSITIVE_INFINITY;
//			int argmin = -1;
//			System.out.println("data = np.array([");
//			for (int j = -20; j <= 20; j ++) {
//				double x = x0[i] + j/10.*Math.abs(dx[i]);
//				xR[i] = x;
//				double y = function.apply(xR);
//				System.out.printf("[%f, %f],\n", x, y);
//				if (y < min) {
//					min = y;
//					argmin = j;
//				}
//			}
//			System.out.printf("]) # %d\n", i);
//			if (Math.abs(argmin) == 20)
//				System.out.println("WARN: check this one; it doesn't seem to have converged.");
		}

		double[][] hessian = new double[x0.length][x0.length]; // then go for the second derivatives
		for (int i = 0; i < hessian.length; i ++) {
			double r = step[i];
			double[] xL = Arrays.copyOf(x0, x0.length);
			xL[i] -= dx[i];
			double l = function.apply(xL);
			if (Double.isFinite(l)) {
				hessian[i][i] = (r - 2*c + l)/(dx[i]*dx[i]); // approximate it as paraboloidic
				for (int j = 0; j < x0.length; j ++) { // and get some diagonal terms
					if (j != i) {
						double u = step[j];
						double[] xUR = Arrays.copyOf(x0, x0.length);
						xUR[i] += dx[i];
						xUR[j] += dx[j];
						double ur = function.apply(xUR);
						hessian[i][j] = hessian[j][i] = (ur - u - r + c)/(dx[i]*dx[j]);
					}
				}
			}
			else { // if we are at a bound
				hessian[i][i] = Math.pow((r - c)/dx[i], 2); // approximate this exponential-ish distribution as gaussian
				for (int j = 0; j < i; j ++)
					hessian[i][j] = hessian[j][i] = 0; // and reset any diagonal terms that previously involved this
			}
		}

		return hessian;
	}

	/**
	 * a discrete representation of an unknown function, capable of evaluating in log time.
	 * 
	 * @author Justin Kunimune
	 */
	public static class DiscreteFunction {
		
		private final boolean equal; // whether the x index is equally spaced
		private final boolean log; // whether to use log interpolation instead of linear
		private final double[] X;
		private final double[] Y;

		/**
		 * instantiate a new function given x and y data in columns. x must
		 * monotonically increase, or the evaluation technique won't work.
		 * @param data array of {x, y}
		 */
		public DiscreteFunction(double[][] data) {
			this(data, false);
		}

		/**
		 * instantiate a new function given x and y data in columns. x must
		 * monotonically increase, or the evaluation technique won't work.
		 * @param data array of {x, y}
		 * @param equal whether the x values are all equally spaced
		 */
		public DiscreteFunction(double[][] data, boolean equal) {
			this(data, equal, false);
		}

		/**
		 * instantiate a new function given x and y data in columns. x must
		 * monotonically increase, or the evaluation technique won't work.
		 * @param data array of {x, y}
		 * @param equal whether the x values are all equally spaced
		 * @param log whether to use log interpolation instead of linear
		 */
		public DiscreteFunction(double[][] data, boolean equal, boolean log) {
			this(collum(data, 0), collum(data, 1));
		}

		/**
		 * instantiate a new function given raw data. x must monotonically
		 * increase, or the evaluation technique won't work.
		 * @param x the x values
		 * @param y the corresponding y values
		 */
		public DiscreteFunction(double[] x, double[] y) {
			this(x, y, false);
		}

		/**
		 * instantiate a new function given raw data. x must monotonically
		 * increase, or the evaluation method won't work.
		 * @param x the x values
		 * @param y the corresponding y values
		 * @param equal whether the x values are all equally spaced
		 */
		public DiscreteFunction(double[] x, double[] y, boolean equal) {
			this(x, y, equal, false);
		}

		/**
		 * instantiate a new function given raw data. x must monotonically
		 * increase, or the evaluation method won't work.
		 * @param x the x values
		 * @param y the corresponding y values
		 * @param equal whether the x values are all equally spaced
		 * @param log whether to use log interpolation instead of linear
		 */
		public DiscreteFunction(double[] x, double[] y, boolean equal, boolean log) {
			if (x.length != y.length)
				throw new IllegalArgumentException("datums lengths must match");
			for (int i = 1; i < x.length; i ++)
				if (x[i] < x[i-1])
					throw new IllegalArgumentException("x must be monotonically increasing.");

			this.X = x;
			this.Y = y;
			this.equal = equal;
			this.log = log;
		}

		/**
		 * it's a function. evaluate it. if this function's x values are equally spaced, this
		 * can be run in O(1) time. otherwise, it will take O(log(n)).
		 * @param x the x value at which to find f
		 * @return f(x)
		 */
		public double evaluate(double x) {
			int i; // we will linearly interpolate x from (X[i], X[i+1]) onto (Y[i], Y[i+1]).
			if (x < X[0]) // if it's out of bounds, we will use the lowest value
				i = 0;
			else if (x >= X[X.length-1]) // or highest value, depending on which is appropriate
				i = X.length-2;
			else if (equal) // nonzero resolution means we can find i itself with linear interpolation
				i = (int)((x - X[0])/(X[X.length-1] - X[0])*(X.length-1)); // linearly interpolate x from X to i
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
			if (log)
				return Y[i]*Math.exp(Math.log(x/X[i])/Math.log(X[i+1]/X[i])*Math.log(Y[i+1]/Y[i]));
			else
				return Y[i] + (x - X[i]) / (X[i+1] - X[i]) * (Y[i+1] - Y[i]); // linearly interpolate x from X[i] to Y[i]
		}
		
		/**
		 * it's a function. evaluate it. if this function's x values are equally spaced, this
		 * can be run in O(1) time. otherwise, it will take O(log(n)).
		 * @param x the x value at which to find f and f's gradient
		 * @return f(x)
		 */
		public Quantity evaluate(Quantity x) {
			int i; // we will linearly interpolate x from (X[i], X[i+1]) onto (Y[i], Y[i+1]).
			if (x.value < X[0]) // if it's out of bounds, we will extrapolate from the lowest values
				i = 0;
			else if (x.value > X[X.length-1]) // or highest values, depending on which is appropriate
				i = X.length-2;
			else if (equal) // nonzero resolution means we can find i itself with linear interpolation
				i = (int)((x.value - X[0])/(X[X.length-1] - X[0])*(X.length-1)); // linearly interpolate x from X to i
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
			if (log)
				return x.over(X[i]).log().times(Math.log(Y[i+1]/Y[i])/Math.log(X[i+1]/X[i])).exp().times(Y[i]);
			else
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
		 * get the second-order-accurate first derivative of this function
		 * @return the derivative as a function of x.
		 */
		public DiscreteFunction derivative() {
			if (!equal)
				throw new IllegalArgumentException("eh");
			int n = X.length;
			double[] dydx = new double[n];
			dydx[0] = (-3*Y[0] + 4*Y[1] - Y[2])/(X[2] - X[0]);
			for (int i = 1; i < n - 1; i ++)
				dydx[i] = (Y[i+1] - Y[i-1])/(X[i+1] - X[i-1]);
			dydx[n-1] = (Y[n-3] - 4*Y[n-2] + 3*Y[n-1])/(X[n-1] - X[n-3]);
			return new DiscreteFunction(X, dydx, true, this.log);
		}
		
		/**
		 * use the trapezoid rule to estimate  the antiderivative of this, with
		 * some arbitrary vertical shift applied.
		 * @return the antiderivative.
		 */
		public DiscreteFunction antiderivative() {
			double[] yOut = new double[X.length];
			yOut[0] = 0; // arbitrarily set the zeroth point to 0
			for (int i = 1; i < X.length; i ++) {
				yOut[i] = yOut[i-1] + (Y[i-1] + Y[i])/2*(X[i] - X[i-1]); // solve for subsequent points using a trapezoid rule
			}
			return new DiscreteFunction(X, yOut, this.equal, this.log);
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

			return new DiscreteFunction(xOut, yOut, true, this.log);
		}

		/**
		 * @return the least x value for which this is not an extrapolation.
		 */
		public double minDatum() {
			return this.X[0];
		}

		/**
		 * @return the greatest value for which this is not an extrapolation.
		 */
		public double maxDatum() {
			return this.X[this.X.length-1];
		}

		@Override
		public String toString() {
			StringBuilder s = new StringBuilder("np.array([");
			for (int i = 0; i < X.length; i ++)
				s.append(String.format(Locale.US, "[%g,%g],", X[i], Y[i]));
			s.append("])");
			return s.toString();
		}
	}
	
	
	/**
	 * A value that tracks its gradient in parameter space for the purpose of error bar
	 * determination.  A unit string may also be attached, but I don't parse or convert
	 * them at all.
	 * @author Justin Kunimune
	 */
	public static class Quantity {
		public final double value;
		public final Vector gradient;
		public final String units;

		public Quantity(double value, int n) {
			this(value, new double[n], "");
		}

		public Quantity(double value, int n, String units) {
			this(value, new double[n], units);
		}

		public Quantity(double value, int i, int n) {
			double[] gradient = new double[n];
			gradient[i] = 1;
			this.value = value;
			this.gradient = new Vector(gradient);
			this.units = "";
		}

		public Quantity(double value, double[] gradient, String units) {
			this(value, new Vector(gradient), units);
		}

		public Quantity(double value, Vector gradient) {
			this(value, gradient, "");
		}

		public Quantity(double value, Vector gradient, String units) {
			this.value = value;
			this.gradient = gradient;
			this.units = units;
		}

		public double variance(double[][] covariance) {
			if (this.getN() == 0)
				return 0;
			if (covariance.length != this.getN() || covariance[0].length != this.getN())
				throw new IllegalArgumentException("this covariance matrix doesn't go with this Quantity");
			
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

		public Quantity exp() {
			return new Quantity(Math.exp(this.value),
			                    this.gradient.times(Math.exp(this.value)));
		}

		public Quantity log() {
			return new Quantity(Math.log(this.value),
			                    this.gradient.times(1/this.value));
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
		 * @return the number of variables upon which this depends
		 */
		public int getN() {
			return this.gradient.getN();
		}

		public Quantity withUnits(String units) {
			return new Quantity(this.value, this.gradient, units);
		}
		
		public String toString(double[][] covariance) {
			if (this.units != null)
				return String.format("%8.6g \u00B1 %8.3g %s",
						this.value, Math.sqrt(this.variance(covariance)), this.units);
			else
				return String.format("%8.6g \u00B1 %8.3g",
						this.value, Math.sqrt(this.variance(covariance)));
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
	
	
	public static void main(String[] args) {
		Random random = new Random();
		for (int i = 0; i < 100; i ++) {
			int n = 2 + (int)(-6*Math.log(Math.random()));
			double p = (Math.cos(Math.random()*Math.PI) + 1)/2.;
			int draw = binomial(n, p, random);
			System.out.printf("%d, %.3f, %d / %.3f = %.3f\n", n, p, draw, n*p, draw/(n*p));
		}
	}
}
