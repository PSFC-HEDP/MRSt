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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.function.BiFunction;
import java.util.function.Function;

/**
 * @author Justin Kunimune
 *
 */
public class Optimization {
	
	/**
	 * perform a backtracking line search with the Armijo-Goldstein condition. the result will
	 * not be a true minimum. it will just be a better point than x0 in the vicinity of x0. it
	 * will not find points beyond x0 ± stepSize.
	 * @param func the function to minimize
	 * @param grad the slope of the function at x0
	 * @param x0 the point to which we are backstepping
	 * @param stepSize the initial step length guess (should be positive)
	 * @return x such that f(x) is better than f(x0)
	 */
	public static double minimizeBacktrack(
			Function<Double, Double> func, double grad, double x0, double stepSize) {
		if (grad > 0)
			stepSize *= -1;
		else if (grad == 0)
			return x0; // if the gradient here is naught, there's absolutely noting we can do
		
		final double c = 0.4, τ = 0.5;
		
		double fx0 = func.apply(x0);
		double x = x0 + stepSize; // take an initial downhill step
		double fx = func.apply(x);
		while (fx - fx0 > (x - x0)*c*grad && Math.abs(x - x0) > Math.abs(stepSize)*1e-15) { // if the function didn't decrease enough (or we hit roundoff error)
			x = x0 + τ*(x - x0); // backstep and try again
			fx = func.apply(x);
		}
		
		return x;
	}
	
	/**
	 * perform a line search with the Wolfe conditions. the result will not be a true minimum,
	 * but it will be near the minimum.
	 * @param func the function to minimize
	 * @param grad the directed gradient function
	 * @param x0 the origin point to which we are backstepping
	 * @param f0 the function at that point (I'm sure you already have it computed)
	 * @param δf0 the slope at that point (I'm not computing it myself)
	 * @param step0 the initial step length guess
	 * @param stepMax the maximum allowed step length
	 * @return x that approximately minimizes f
	 * @throws RuntimeException if the function is not smooth or otherwise ABNO
	 */
	public static double minimizeWolfe(
			Function<Double, Double> func, Function<Double, Double> grad,
			double x0, double f0, double δf0, double step0, double stepMax) throws RuntimeException {
		if (step0 <= 0)
			throw new IllegalArgumentException("Initial step must be positive.");
		if (step0 > stepMax)
			throw new IllegalArgumentException("Initial step must be bounded by stepMax.");
		if (δf0 > 0)
			throw new IllegalArgumentException("Initial step must be downhill.");
		if (δf0 == 0)
			return x0; // if the gradient here is naught, there's absolutely noting we can do
		if (!Double.isFinite(f0))
			throw new IllegalArgumentException("Initial guess was "+f0);
		
		final double α = 1e-4, β = 0.9;
		
		double med = x0 + step0; // take an initial downhill step
		double min = x0, max = x0 + 2*stepMax;
		double truMax = x0 + stepMax;
		double lowestPlace = x0, lowestValue = f0;
		if (Double.isNaN(lowestValue))
			throw new IllegalStateException("You little bitch.");
		while (true) {
			if (Double.isNaN(lowestValue))
				throw new IllegalStateException("How?! When?!");
			if ((max - min)/Math.max(med, step0) < 1e-10) { // if this has become quite tight
				if (Double.isNaN(func.apply(lowestPlace)))
					throw new IllegalStateException("waaaa how did this happen "+lowestPlace+" and "+lowestValue+" and "+med);
				return lowestPlace;
			}
			med = Math.min(med, truMax); // enforce that it not go past its true maximum
			double f = func.apply(med);
			if (Double.isNaN(f)) {
				max = med; // ew, get away from it!
				med = (max + min)/2;
				continue;
			}
			if (f < lowestValue) { // keep track of the lowest value we could find just in case all fails
				lowestPlace = med;
				lowestValue = f;
				if (Double.isNaN(lowestValue))
					throw new IllegalStateException("I literally just asked you if f was NaN");
			}
			if (f > f0 + α*(med - x0)*δf0) { // if the decrease condition is not met
				max = med; // we need to go closer
				med = (max + min)/2;
				continue;
			}
			double δf = grad.apply(med);
			if (Math.abs(δf) > β*Math.abs(δf0)) { // if the curvature condition is not met
				if (δf > 0) {
					max = med; // we might need to go closer
					med = (max + min)/2;
				}
				else if (med == truMax) {
//					if (Double.isNaN(func.apply(truMax)))
//						throw new IllegalStateException("But I already asked if f(med) was NaN! "+f);
					return truMax; // this might be the best we can get
				}
				else {
					min = med; // but most likely the local min is farther out
					med = Math.min((max + min)/2, 2.7*min);
				}
				continue;
			}
			
//			if (Double.isNaN(func.apply(med)))
//				throw new IllegalArgumentException("What are you doing? I asked you if f was NaN and you sed no!");
			return med; // if both are met, we're done here
		}
	}
	
	
	/**
	 * perform Nelder-Mead unconstrained optimization, making a SWAG for the initial simplex
	 * scale.
	 * @param func the function to minimize
	 * @param x0 the initial guess of the minimizing parameter vector
	 * @param tol the relative tolerance for error
	 * @return the x vector that minimizes f(x)
	 */
	public static double[] minimizeNelderMead(
			Function<double[], Double> func, double[] x0, double tol) {
		double[] scale = new double[x0.length];
		for (int i = 0; i < x0.length; i ++)
			scale[i] = Math.max(1, Math.abs(x0[i])); // guess the scale based on the initial guess
		return minimizeNelderMead(func, x0, scale, tol);
	}
	
	
	/**
	 * perform Nelder-Mead unconstrained optimization.
	 * @param func the function to minimize
	 * @param x0 the initial guess of the minimizing parameter vector
	 * @param scale an array of the same size as x0 giving the general length scale of its
	 * derivatives, to kick off the simplex.
	 * @param tol the relative tolerance for error
	 * @return the x vector that minimizes f(x)
	 */
	public static double[] minimizeNelderMead(
			Function<double[], Double> func, double[] x0, double[] scale, double tol) {
		final double α = 1, γ = 2, ρ = 1/2., σ = 1/2.; // declare values
		
		Function<Matrix, Double> f = (mat) -> func.apply(mat.values[0]); // change these from double[] to Matrix things
		
		Matrix[] x = new Matrix[x0.length+1]; // the vertices of the simplex
		double[] fx = new double[x0.length+1]; // the values of f at the vertices
		for (int i = 0; i < x.length; i ++) {
			x[i] = new Matrix(new double[][] {x0.clone()}); // initialize the vertices as the guess
			if (i < x0.length)
				x[i].set(0, i, x[i].get(0, i) + scale[i]/6); // with multidimensional perturbations
			fx[i] = f.apply(x[i]); // and get the initial values
		}
		if (!Double.isFinite(fx[x0.length]))
			throw new IllegalArgumentException("Initial guess yielded bunk value: "+fx[x0.length]);
		
		while (true) { // now for the iterative part
			int iWorst = NumericalMethods.argmax(fx);
			int iNext = NumericalMethods.argpenmax(fx);
			int iBest = NumericalMethods.argmin(fx);
			
			boolean done = true;
			for (int j = 0; j < scale.length; j ++)
				if (Math.abs(x[iBest].get(0,j) - x[iWorst].get(0,j)) > tol*scale[j])
					done = false;
			if (done)
				return x[iBest].values[0];
			
			Matrix xC = new Matrix(1, x0.length);
			for (int i = 0; i < x.length; i ++) // compute the best-guess centroid
				if (i != iWorst)
					xC = xC.plus(x[i].over(x.length - 1));
			
			Matrix xR = xC.plus(xC.minus(x[iWorst]).times(α)); // compute the reflected point
			double fxR = f.apply(xR);
			
			if (fxR < fx[iBest]) { // if this is the best point yet
				Matrix xE = xC.plus(xR.minus(xC).times(γ)); // compute the expanded point
				double fxE = f.apply(xE);
				
				if (fxE < fxR) {
					x[iWorst] = xE;
					fx[iWorst] = fxE;
				}
				else {
					x[iWorst] = xR;
					fx[iWorst] = fxR;
				}
			}
			else if (fxR < fx[iNext]) { // if this is better than the second worst
				x[iWorst] = xR;
				fx[iWorst] = fxR;
			}
			else if (fxR < fx[iWorst]) { // if this has just a little bit of redeeming quality
				Matrix xS = xC.plus(xR.minus(xC).times(ρ));
				double fxS = f.apply(xS);
				
				if (fxS <= fxR) {
					x[iWorst] = xS;
					fx[iWorst] = fxS;
				}
				else {
					for (int i = 0; i < x.length; i ++) {
						if (i != iBest) {
							x[i] = x[iBest].plus(x[i].minus(x[iBest]).times(σ)); // move all vertices inward
							fx[i] = f.apply(x[i]);
						}
					}
				}
			}
			else {
				Matrix xS = xC.plus(x[iWorst].minus(xC).times(ρ)); // compute the contracted point
				double fxS = f.apply(xS);
				
				if (fxS < fx[iWorst]) {
					x[iWorst] = xS;
					fx[iWorst] = fxS;
				}
				else {
					for (int i = 0; i < x.length; i ++) {
						if (i != iBest) {
							x[i] = x[iBest].plus(x[i].minus(x[iBest]).times(σ)); // move all vertices inward
							fx[i] = f.apply(x[i]);
						}
					}
				}
			}
		}
	}
	
	
	/**
	 * perform a simple coordinate descent optimization scheme without a provided gradient.
	 * practically, the gradient will be estimated by evaluating the function with
	 * perturbations.
	 * @param func the function to be optimized
	 * @param x0 the initial guess
	 * @param tol the relative tolerance for error
	 * @return the x that minimizes func
	 */
	public static double[] minimizeCoordinateDescent(
			Function<double[], Double> func, double[] x0, double tol) {
		double[] scale = new double[x0.length];
		for (int i = 0; i < x0.length; i ++)
			scale[i] = Math.max(1, Math.abs(x0[i]))*.7; // guess the scale based on the initial guess
		return minimizeCoordinateDescent(func, x0, scale, tol);
	}
	
	
	/**
	 * perform a simple coordinate descent optimization scheme without a provided gradient.
	 * practically, the gradient will be estimated by evaluating the function with
	 * perturbations.
	 * @param func the function to be optimized
	 * @param x0 the initial guess
	 * @param scale an array of the same size as x0 giving the general length scale of its
	 * derivatives, to kick off the line search.
	 * @param tol the relative tolerance for error
	 * @return the x that minimizes func
	 */
	public static double[] minimizeCoordinateDescent(
			Function<double[], Double> func, double[] x0, double[] scale, double tol) {
		final double dt = 1e-9;
		List<Function<double[], Double>> grad = new ArrayList<Function<double[], Double>>(x0.length);
		for (int i = 0; i < x0.length; i ++) {
			final int I = i;
			grad.add((x) -> {
				double y0 = func.apply(x);
				x[I] += scale[I]*dt;
				double y1 = func.apply(x);
				x[I] -= scale[I]*dt;
				return (y1 - y0)/dt;
			});
		}
		return minimizeCoordinateDescent(func, grad, x0, scale, tol);
	}
	
	
	/**
	 * perform a simple coordinate descent optimization scheme.
	 * @param func the function to be optimized
	 * @param grad an array of derivatives, each corresponding to one dimension: df(x)/dxi
	 * @param x0 the initial guess
	 * @param scale an array of the same size as x0 giving the general length scale of its
	 * derivatives, to kick off the line search.
	 * @param tol the relative tolerance for error
	 * @return the x that minimizes func
	 */
	public static double[] minimizeCoordinateDescent(
			Function<double[], Double> func, List<Function<double[], Double>> grad,
			double[] x0, double[] scale, double tol) {
		double[] x = x0.clone();
		double fx = func.apply(x);
		if (!Double.isFinite(fx))
			throw new IllegalArgumentException("Initial guess yielded bunk value");
		
		while (true) {
			for (int i = 0; i < x0.length; i ++) {
				final int I = i;
				x[i] = minimizeBacktrack((xi) -> {
					double[] xp = x.clone();
					xp[I] = xi;
					return func.apply(xp);
				}, grad.get(i).apply(x), x[i], scale[i]);
				System.out.println(-func.apply(x));
			}
			double fxp = func.apply(x);
			if (Math.abs((fxp - fx)/fx) < tol)
				return x;
			else
				fx = fxp;
		}
	}
	
	
	/**
	 * perform a coordinate descent Nelder Mead hybrid. it treats clusters of n coordinates as
	 * single dimensions, and optimizes along each of those macrodimensions individually using
	 * Nelder-Mead. good for cases where coordinates cluster into relatively unrelated sets
	 * (i.e. the Hessian is blocky or something). idk if this is actually useful for anything.
	 * @param func the function to minimize
	 * @param x0 the initial guess
	 * @param n the number of dimensions per macrodimension
	 * @param tol the relative tolerance for error
	 * @return the value of x that minimizes func
	 */
	public static double[] minimizeCoordinateDescent(
			Function<double[], Double> func, double[] x0, int n, double tol) {
		double[] x = x0.clone();
		double fx = func.apply(x);
		while (true) {
			for (int i = 0; i < x.length; i += n) { // iterate through the coordinates
				double[] y0 = Arrays.copyOfRange(x, i, Math.min(i+n, x.length)); // extract slices of the array
				final int I = i;
				double[] yOpt = minimizeNelderMead((y) -> { // then run Nelder-Mead on that slice
					double[] xp = x.clone();
					System.arraycopy(y, 0, xp, I, y.length);
					return func.apply(xp);
				}, y0, tol);
				System.arraycopy(yOpt, 0, x, i, yOpt.length);
				System.out.println(Arrays.toString(x)+",");
			} // after going through each coordinate
			double fxp = func.apply(x); // check how much change that whole cycle got us
			if (Math.abs((fxp - fx)/fx) < tol)
				return x;
			else
				fx = fxp;
		}
	}
	
	
	/**
	 * perform L-BFGS-B constrained optimization without a provided gradient. practically, the
	 * gradient will be estimated with finite differences.
	 * <p>
	 * Richard H. Byrd, Peihuang Lu, Jorge Nocedal, & Ciyou Zhu (1994). "A limited memory
	 *   algorithm for bound constrained optimization." <i>Northwestern University Department
	 *   of Electrical Engineering and Computer Science.</i> Technical Report NAM-08.
	 * @param func the function to minimize
	 * @param x0 the initial guess of the minimizing parameter vector
	 * @param scale the approximate scale at which the funccion should vary with each variable
	 * @param lower the lowest allowable values of each of the elements of x
	 * @param upper the greatest allowable values of each of the elements of x
	 * @param relTol the relative tolerance for error (set to 0 to use only absolute)
	 * @param absTol the absolute tolerance for error (set to 0 to use only relative)
	 * @return the vector x that minimizes f(x)
	 */
	public static double[] minimizeLBFGSB(
			Function<double[], Double> func, double[] x0, double[] scale, double[] lower, double[] upper,
			double relTol, double absTol) {
		if (x0.length != scale.length)
			throw new IllegalArgumentException("the scale must be the same size as the other dimensional things!");
		for (double value: scale)
			if (value == 0 || Double.isInfinite(value))
				throw new IllegalArgumentException("all scales must be finite scalars, not " + value);

		final double ds = 1e-5;
		Function<double[], double[]> gradient = (x) -> { // finite difference gradient:
			double y0 = func.apply(x);
			double[] dydx = new double[x.length];
			for (int i = 0; i < x.length; i ++) { // compute each dimension individually
				double dx = scale[i]*ds;
				if (x[i] + dx >= upper[i]) dx = -dx; // chose a direction to avoid the bounds
				x[i] += dx; // perturb
				dydx[i] = (func.apply(x) - y0)/dx; // measure
				x[i] -= dx; // unperturb
			}
			return dydx;
		};
		BiFunction<double[], double[], Double> linGradient = (x, v) -> { // finite difference linear gradient:
			double y0 = func.apply(x);
			double dx = ds; // the step direction must have the same sign for all components now
			for (int i = 0; i < x.length; i ++)
				if (x[i] + v[i]*dx > upper[i] || x[i] + v[i]*dx < lower[i])
					dx = -dx; // but we can still flip it if one of the components is near a bound
			for (int i = 0; i < x.length; i ++)
				x[i] += v[i]*dx; // do all the dimensionless perturbations
			double dydx = (func.apply(x) - y0)/dx; // measure
			for (int i = 0; i < x.length; i ++)
				x[i] -= v[i]*dx; // then put them all back
			return dydx; // the result is v-dot-grady
		};
		return minimizeLBFGSB(func, gradient, linGradient, x0, lower, upper, relTol, absTol);
	}
	
	
	/**
	 * perform L-BFGS-B constrained optimization with a provided gradient.
	 * <p>
	 * Richard H. Byrd, Peihuang Lu, Jorge Nocedal, & Ciyou Zhu (1994). "A limited memory
	 *   algorithm for bound constrained optimization." <i>Northwestern University Department
	 *   of Electrical Engineering and Computer Science.</i> Technical Report NAM-08.
	 * @param func the function to minimize
	 * @param grad the gradient of that function
	 * @param x0 the initial guess of the minimizing parameter vector
	 * @param lower the lowest allowable values of each of the elements of x
	 * @param upper the greatest allowable values of each of the elements of x
	 * @param relTol the relative tolerance for error (set to 0 to use only absolute)
	 * @param absTol the absolute tolerance for error (set to 0 to use only relative)
	 * @return the vector x that minimizes f(x)
	 */
	public static double[] minimizeLBFGSB(
			Function<double[], Double> func, Function<double[], double[]> grad, double[] x0,
			double[] lower, double[] upper, double relTol, double absTol) {
		return minimizeLBFGSB(
				func, grad, (x, v) -> {
					double[] g = grad.apply(x);
					double gv = 0;
					for (int i = 0; i < v.length; i ++)
						gv += g[i]*v[i];
					return gv;
				}, x0, lower, upper, relTol, absTol);
	}
	
	/**
	 * Perform L-BFGS-B constrained optimization with a provided multidimensional gradient. In
	 * addition, since it may be much faster to compute the gradient in a single direction,
	 * this function also accepts a separate linear gradient function.
	 * <p>
	 * Richard H. Byrd, Peihuang Lu, Jorge Nocedal, & Ciyou Zhu (1994). "A limited memory
	 *   algorithm for bound constrained optimization." <i>Northwestern University Department
	 *   of Electrical Engineering and Computer Science.</i> Technical Report NAM-08.
	 * @param func the function to minimize
	 * @param grad the gradient of that function
	 * @param linGrad the gradient of that function dotted with the second input
	 * @param x0 the initial guess of the minimizing parameter vector
	 * @param lower the lowest allowable values of each of the elements of x
	 * @param upper the greatest allowable values of each of the elements of x
	 * @param relTol the relative tolerance for error (set to 0 to use only absolute)
	 * @param absTol the absolute tolerance for error (set to 0 to use only relative)
	 * @return the vector x that minimizes f(x)
	 */
	public static double[] minimizeLBFGSB(
			Function<double[], Double> func, Function<double[], double[]> grad, BiFunction<double[], double[], Double> linGrad,
			double[] x0, double[] lower, double[] upper, double relTol, double absTol) {
		if (lower.length != x0.length || upper.length != x0.length)
			throw new IllegalArgumentException("Bound lengths don't match");
		for (int i = 0; i < x0.length; i ++)
			if (lower[i] > x0[i] || x0[i] > upper[i])
				throw new IllegalArgumentException("Upper bounds must be greater than lower bounds");
		for (double value: x0)
			if (Double.isNaN(value))
				throw new IllegalArgumentException("NaNs are strictly forbidden.");
		
		final int mMax = 6;
		final int n = x0.length;
		
		Function<Matrix, Double> funcMat = (mat) -> func.apply(mat.T().values[0]); // change these from double[] to Matrix things
		Function<Matrix, Matrix> gradMat = (mat) -> new Matrix(new double[][] {grad.apply(mat.T().values[0])}).T();
		BiFunction<Matrix, Matrix, Double> linGradMat = (mat, v) -> linGrad.apply(mat.T().values[0], v.T().values[0]);
		
		int iter = 0;
		LinkedList<Matrix> yHist = new LinkedList<Matrix>();
		LinkedList<Matrix> sHist = new LinkedList<Matrix>();
		double θ = 1;
		Matrix xk = new Matrix(new double[][] {x0.clone()}).T(); // the current best guess (column matrix)
		Matrix gk = gradMat.apply(xk); // the current value of g from the L-BFGS algorithm
		double fxk = funcMat.apply(xk);
		if (!Double.isFinite(fxk))
			throw new IllegalArgumentException("Initial guess yielded bunk value");
		
		while (true) {
//			if (Math.random() < 1e-2)
//				System.out.println("FINER: "+fxk);
			
			final int m = yHist.size();
			Matrix dk;
			if (m > 0) {
				assert θ > 0;
				for (int i = 0; i < m; i ++)
					assert sHist.get(i).dot(yHist.get(i)) > 0;
				Matrix Dk = new Matrix(m, m); // STEP 1: approximate the inverse Hessian
				for (int i = 0; i < m; i ++)
					Dk.set(i, i, yHist.get(i).dot(sHist.get(i)));
				Matrix Wk = new Matrix(n, 2*m);
				for (int i = 0; i < n; i ++) {
					for (int j = 0; j < m; j ++) {
						Wk.set(i, j,   yHist.get(j).get(i, 0));
						Wk.set(i, j+m, θ*sHist.get(j).get(i, 0));
					}
				}
				Matrix Lk = new Matrix(m, m);
				for (int i = 0; i < m; i ++)
					for (int j = 0; j < i; j ++)
						Lk.set(i, j, sHist.get(i).dot(yHist.get(j)));
				Matrix Mkinv = new Matrix(2*m, 2*m);
				for (int i = 0; i < m; i ++) {
					Mkinv.set(i, i, -Dk.get(i, i));
					for (int j = 0; j < m; j ++) {
						Mkinv.set(i+m, j, Lk.get(i, j));
						Mkinv.set(i, j+m, Lk.get(j, i));
						Mkinv.set(i+m, j+m, θ*sHist.get(i).dot(sHist.get(j)));
					}
				}
				Matrix Mk = Mkinv.inv();
				for (int i = 0; i < Mk.getN(); i ++)
					for (int j = 0; j < Mk.getM(); j ++)
						if (Double.isNaN(Mk.get(i, j))) {
							System.out.println(θ);
							System.out.println(yHist.getLast().dot(yHist.getLast())/yHist.getLast().dot(sHist.getLast()));
							System.out.println(yHist.getLast().dot(sHist.getLast())/sHist.getLast().dot(sHist.getLast()));
							System.out.println(Mkinv);
							System.out.println(Mk);
							for (int ih = 0; ih < Mk.getN(); ih ++) {
								for (int jo = 0; jo < Mk.getM(); jo ++) {
									if (ih != jo)
										Mkinv.set(ih, jo, Mkinv.get(ih, jo)/2);
								}
							}
							Mk = Mkinv.inv();
						}

				double[] breakpoints = new double[n]; // STEP 2: find the Cauchy point -- the quadratic minimum in the downhill direction
				Matrix d = new Matrix(n, 1);
				for (int i = 0; i < n; i ++) {
					if (gk.get(i, 0) < 0)
						breakpoints[i] = (xk.get(i, 0) - upper[i])/gk.get(i, 0);
					else if (gk.get(i, 0) > 0)
						breakpoints[i] = (xk.get(i, 0) - lower[i])/gk.get(i, 0);
					else
						breakpoints[i] = Double.POSITIVE_INFINITY;
					assert breakpoints[i] >= 0;
					d.set(i, 0, (breakpoints[i] == 0) ? 0 : -gk.get(i, 0));
				}
				if (d.norm() == 0) // if there is no step here
					return xk.T().values[0]; // there's nothing to do; just return what you have
				
				LinkedList<Integer> breakpointOrder = new LinkedList<Integer>();
				for (int i = 0; i < n; i ++) // in terms of the original paper, this list contains the order in which indices are removed from F
					if (breakpoints[i] > 0)
						breakpointOrder.add(i);
				breakpointOrder.sort((iA, iB) -> (int)Math.signum(breakpoints[iA] - breakpoints[iB]));
				
				Matrix p = Wk.T().times(d); // first look for minimum in first segment
				Matrix c = new Matrix(2*m, 1);
				double told = 0;
				int b = breakpointOrder.pop(); // the index that hits bound at the end of this interval
				double t = breakpoints[b];
				double Δt = t;
				double dfdt = d.dot(gk);
				double d2fdt2 = -θ*dfdt - p.dot(Mk.times(p));
				assert dfdt < 0 : d;
				double Δtmin = (d2fdt2 != 0) ? -dfdt/d2fdt2 : Δt;
				while (Double.isFinite(Δt) && Δtmin >= Δt) { // then check all subsequent segments
					if (breakpointOrder.size() == 0) {
						Δtmin = Δt;
						break;
					}
					double xCb = (d.get(b, 0) > 0) ? upper[b] : lower[b];
					double zb = xCb - xk.get(b, 0);
					double gb = gk.get(b, 0);
					Matrix wb = Wk.getRow(b);
					c = c.plus(p.times(Δt));
					d.set(b, 0, 0);
					dfdt = dfdt + Δt*d2fdt2 + gb*gb + θ*gb*zb - gb*wb.dot(Mk.times(c));
					d2fdt2 = d2fdt2 - θ*gb*gb - 2*gb*wb.dot(Mk.times(p)) - gb*gb*wb.dot(Mk.times(wb));
					if (d2fdt2 < 0) {
						System.err.println("WARN: encountered a concave-down line search");
						d2fdt2 = 0;
					}
					p = p.plus(wb.times(gb));
					told = t;
					b = breakpointOrder.pop();
					t = breakpoints[b];
					Δt = t - told;
					if (d2fdt2 != 0)
						Δtmin = -dfdt/d2fdt2;
					else if (dfdt > 0)
						Δtmin = 0;
					else
						Δtmin = Δt;
					assert !Double.isNaN(Δtmin): Δtmin+"\n"+dfdt+"\n"+d2fdt2+"\n"+Δt;
				}
				if (Δtmin < 0)
					Δtmin = 0;
				else if (Double.isInfinite(Δtmin))
					Δtmin = 1;
				double tC = told + Δtmin;
				assert Double.isFinite(tC): tC+"\n"+told+"\n"+Δtmin+"\n"+Arrays.toString(breakpoints);
				c = c.plus(p.times(Δtmin));
				List<Integer> F = new ArrayList<Integer>(breakpointOrder.size()+1);
				Matrix xC = new Matrix(n, 1);
				for (int i = 0; i < n; i ++) { // I'm setting xC here not how it's done in the paper, because the paper version is _totally_ wrong for this part
					if (breakpoints[i] <= tC)
						xC.set(i, 0, (gk.get(i, 0) < 0) ? upper[i] : lower[i]);
					else
						xC.set(i, 0, xk.get(i, 0) - tC*gk.get(i, 0));
					if (xC.get(i, 0) != lower[i] && xC.get(i, 0) != upper[i])
						F.add(i);
				}
				assert xC.equals(P(xC, lower, upper)) : tC+",\n"+xk+",\n"+gk+"\n"+xC+",\n"+P(xC, lower, upper);
				
				if (F.size() >= 1) {
					Matrix Zk = new Matrix(n, F.size()); // STEP 3: find an approximate bound minimum (direct primal method)
					for (int f = 0; f < F.size(); f ++)
						Zk.set(F.get(f), f, 1);
					Matrix rHatC = Zk.T().times(gk.plus(xC.minus(xk).times(θ)).minus(Wk.times(Mk.times(c))));
					Matrix v = Wk.T().times(Zk.times(rHatC));
					v = Mk.times(v);
					Matrix N = Mk.times(Wk.T().times(Zk.times(Zk.T().times(Wk)))).over(-θ);
					for (int i = 0; i < 2*m; i ++)
						N.set(i, i, 1 + N.get(i, i));
					v = N.inv().times(v);
					Matrix dHatU = rHatC.over(θ).plus(Zk.T().times(Wk.times(v)).over(θ*θ)).times(-1); // the paper has a sign error
//					assert dHatU.dot(rHatC) < 0;
					if (dHatU.dot(rHatC) > 0) // this part is tricky; make sure you're stepping downhill from the Cauchy point
						dHatU = dHatU.times(-1);
					Matrix dU = Zk.times(dHatU);
					double αStar = 1;
					for (int i = 0; i < n; i ++) {
						if (xC.get(i, 0) + dU.get(i, 0) > upper[i]) // the paper has an error here; you need to multiply dU by Zk before comparing to u and l
							αStar = Math.min(αStar, (upper[i] - xC.get(i, 0))/dU.get(i, 0));
						else if (xC.get(i, 0) + dU.get(i, 0) < lower[i])
							αStar = Math.min(αStar, (lower[i] - xC.get(i, 0))/dU.get(i, 0));
					}
					assert αStar > 0 && αStar <= 1 : αStar;
					Matrix xBar = xC.plus(dU.times(αStar));
					xBar = P(xBar, lower, upper); // strictly speaking I shouldn't need this, but roundoff
					dk = xBar.minus(xk);
				}
				else {
					dk = xC.minus(xk);
				}
			}
			else {
				dk = P(xk.minus(gk), lower, upper).minus(xk);
			}
			if (gk.dot(dk) > 0) {
				System.err.println("WARN: it tried to step uphill.");
				dk = P(xk.minus(dk), lower, upper).minus(xk);
			}
			if (gk.dot(dk) > 0) {
				System.err.println("WARN: wait, it's all uphill?");
				return P(xk, lower, upper).T().values[0]; // well, dang.
			}
			
			double λMax = Double.POSITIVE_INFINITY; // STEP 4: perform a line search
			for (int i = 0; i < n; i ++) { // enforcing λMax to ensure the enhanced steps stay in bounds
				double λBi = Double.POSITIVE_INFINITY;
				if (dk.get(i, 0) > 0)
					λBi = (upper[i] - xk.get(i, 0))/dk.get(i, 0);
				else if (dk.get(i, 0) < 0)
					λBi = (lower[i] - xk.get(i, 0))/dk.get(i, 0);
				if (λBi < λMax)
					λMax = λBi;
			}
			assert λMax >= 1 : "Why is lambda max imposing "+λMax;
			final Matrix Xk = xk, Dk = dk;
			double λk = minimizeWolfe(
					(λ) -> {
						return funcMat.apply(P(Xk.plus(Dk.times(λ)), lower, upper)); // strictly speaking I shouldn't need to call P here, but there's roundoff stuff
					}, (λ) -> {
						return linGradMat.apply(P(Xk.plus(Dk.times(λ)), lower, upper), Dk);
					}, 0, fxk, gk.dot(dk), 1, λMax);
			xk = P(xk.plus(dk.times(λk)), lower, upper);
			for (int i = 0; i < xk.getN(); i ++) {
				if (Double.isNaN(xk.get(i, 0))) {
					throw new IllegalArgumentException("whence did this NaN come from??");
				}
			}

			Matrix gkp1 = gradMat.apply(xk);
			double fxkp1 = funcMat.apply(xk);
			if (Double.isNaN(fxkp1))
				throw new IllegalArgumentException("no, no. You can't do that.");
			if ((fxk - fxkp1)/Math.abs(fxk) <= relTol || fxk - fxkp1 <= absTol) { // STEP 5: stop condition
				return xk.T().values[0]; // if we're into it and the energy isn't really changing, then we're done
			}
			else if (iter > 2000) {
				System.err.println("WARN: Maximum iterations reached.");
				return xk.T().values[0];
			}
			assert fxkp1 <= fxk: fxkp1+", "+fxk+", "+Math.abs(fxk)+", "+relTol+", "+absTol;
			
			Matrix sk = dk.times(λk); // STEP 6: save historical vector information
			Matrix yk = gkp1.minus(gk);
			if (yk.dot(sk) > 1e-15*yk.norm()) {
				yHist.addLast(yk);
				sHist.addLast(sk);
				θ = yk.dot(sk)/sk.dot(sk);
			}
			if (yHist.size() > mMax) {
				yHist.removeFirst();
				sHist.removeFirst();
			}
			gk = gkp1;
			fxk = fxkp1;
			iter ++;
		}
	}
	
	private static Matrix P(Matrix x, double[] l, double[] u) {
		Matrix y = new Matrix(x.getN(), 1);
		for (int i = 0; i < x.getN(); i ++)
			y.set(i, 0, Math.max(Math.min(x.get(i, 0), u[i]), l[i]));
		return y;
	}
	
	
	/**
	 * perform L-BFGS unconstrained optimization without a provided gradient. practically, the
	 * gradient will be estimated by evaluating the function with perturbations.
	 * @param func the function to minimize
	 * @param x0 the initial guess of the minimizing parameter vector
	 * @param tol the relative tolerance for error
	 * @return the x vector that minimizes f(x)
	 */
	public static double[] minimizeLBFGS(
			Function<double[], Double> func, double[] x0, double tol) {
		final double dx = 1e-9;
		Function<double[], double[]> gradient = (x) -> {
			double y0 = func.apply(x);
			double[] dydx = new double[x.length];
			for (int i = 0; i < x.length; i ++) {
				x[i] += dx;
				dydx[i] = (func.apply(x) - y0)/dx;
				x[i] -= dx;
			}
			return dydx;
		};
		return minimizeLBFGS(func, gradient, x0, tol);
	}
	
	
	/**
	 * perform L-BFGS unconstrained optimization.
	 * @param func the function to minimize
	 * @param grad the gradient of the function
	 * @param x0 the initial guess of the minimizing parameter vector
	 * @param tol the relative tolerance for error
	 * @return the x vector that minimizes f(x)
	 */
	public static double[] minimizeLBFGS(
			Function<double[], Double> func, Function<double[], double[]> grad, double[] x0, double tol) {
		final double M = 6;
		
		Function<Matrix, Double> funcMat = (mat) -> func.apply(mat.T().values[0]); // change these from double[] to Matrix things
		Function<Matrix, Matrix> gradMat = (mat) -> new Matrix(new double[][] {grad.apply(mat.T().values[0])}).T();
		
		LinkedList<Matrix> sHist = new LinkedList<Matrix>();
		LinkedList<Matrix> yHist = new LinkedList<Matrix>();
		Matrix gkMinus1 = null; // the previous value of $g$ from the L-BFGS algorithm
		Matrix x = new Matrix(new double[][] {x0.clone()}).T(); // the current best guess (column matrix)
		
		double Ui = funcMat.apply(x);
		if (!Double.isFinite(Ui))
			throw new IllegalArgumentException("Initial guess yielded bunk value");
		
		while (true) {
			Matrix gk = gradMat.apply(x);
			
			if (gkMinus1 != null) // STEP 5 (cont.): save historical vector information
				yHist.addLast(gk.minus(gkMinus1));
			if (sHist.size() > M) {
				sHist.removeFirst();
				yHist.removeFirst();
			}
			
			Matrix q = gk.times(-1); // STEP 2: choose the step direction
			double[] alpha = new double[sHist.size()];
			for (int i = sHist.size()-1; i >= 0; i --) { // this is where it gets complicated
				alpha[i] = sHist.get(i).dot(q)/yHist.get(i).dot(sHist.get(i)); // see the paper cited at the top, page 779.
				q = q.plus(yHist.get(i).times(-alpha[i]));
			}
			double H0;
			if (!sHist.isEmpty())
				H0 = sHist.getLast().dot(yHist.getLast())/yHist.getLast().dot(yHist.getLast()); // this is our very rough estimate of the inverse Hessian
			else
				H0 = 1/q.norm();
			Matrix dk = q.times(H0);
			for (int i = 0; i < sHist.size(); i ++) {
				double beta = yHist.get(i).dot(dk)/yHist.get(i).dot(sHist.get(i));
				dk = dk.plus(sHist.get(i).times(alpha[i]-beta));
			}
			
			double gradDotVel = gk.dot(dk);
			if (gradDotVel > 0) { // ensure this number is never positive
				System.err.printf("WARN: It tried to step uphill with g_k \\cdot d_k = %f. I don't know what that means.\n", gradDotVel);
				dk = dk.times(-1);
				gradDotVel = gk.dot(dk);
			}
			
			final Matrix X = x, Dk = dk;
			double timestep = minimizeBacktrack((dt) -> { // STEP 3: choose the step size
				return funcMat.apply(X.plus(Dk.times(dt)));
			}, gradDotVel, 0, 1);
			x = x.plus(dk.times(timestep));
			double Uf = funcMat.apply(x);
			
			if (sHist.size() == M && Math.abs((Ui - Uf)/Ui) < tol) { // STEP 4: stop condition
				return x.T().values[0]; // if we're into it and the energy isn't really changing, then we're done
			}
			
			Matrix sk = dk.times(timestep); // STEP 5: save historical vector information
			sHist.addLast(sk);
			gkMinus1 = gk;
			Ui = Uf;
		}
	}
	
	
	/**
	 * A special case of Gelfgat et al.'s deconvolution scheme. Images are rectangular and the
	 * input and output have the same size.
	 * @param F n×m desired image
	 * @param D n×m measurement variance matrix
	 * @param P nm×nm transfer matrix
	 * @param alpha the smoothing
	 */
	public static double[][] optimizeGelfgat(double[][] F, double[][] D, double[][] P, double alpha) {
		if (F.length != D.length || F[0].length != D[0].length || P.length != P[0].length || P.length != D.length*D[0].length)
			throw new IllegalArgumentException("I can't work with this; have you seen these dimensions‽ "+F.length+"×"+F[0].length+", "+D.length+"×"+D[0].length+", "+P.length+"×"+P[0].length+"!");
		final int n = F.length, m = F[0].length;
		
		double[][] g = new double[n][m];
		for (int i = 0; i < n; i ++)
			for (int j = 0; j < m; j ++)
				g[i][j] = 1.;
		double G;
		
		int iter = 0;
		double score = Double.NEGATIVE_INFINITY, scorePrev;
		do { // use Gelfgat et al.'s program to deconvolve the spectrum
			double Σg = 0;
			for (int i = 0; i < n; i ++)
				for (int j = 0; j < m; j ++)
					Σg += g[i][j];
			for (int i = 0; i < n; i ++)
				for (int j = 0; j < m; j ++)
					g[i][j] /= Σg; // renormalize g to account for roundoff
			
			double[][] s = gelfgatConvolve(P, g);
			
			double ΣFs = 0, Σss = 0;
			for (int k = 0; k < n; k ++) {
				for (int l = 0; l < m; l ++) {
					ΣFs += F[k][l]*s[k][l]/D[k][l];
					Σss += s[k][l]*s[k][l]/D[k][l];
				}
			}
			G = ΣFs/Σss; // compute the optimal G
			
			double[][] δg = new double[n][m];
			for (int i = 0; i < n; i ++)
				for (int j = 0; j < m; j ++)
					for (int k = 0; k < n; k ++)
						for (int l = 0; l < m; l ++)
							δg[i][j] += g[i][j] * P[m*k+l][m*i+j]*(F[k][l] - G*s[k][l])/D[k][l];
			
			double[][] δs = gelfgatConvolve(P, δg);
			
			double Fδ = 0, Ss = 0, Sδ = 0, Dδ = 0;
			for (int k = 0; k < n; k ++) {
				for (int l = 0; l < m; l ++) {
					Fδ += F[k][l]*δs[k][l]/D[k][l];
					Ss += s[k][l]*s[k][l]/D[k][l];
					Sδ += s[k][l]*δs[k][l]/D[k][l];
					Dδ += δs[k][l]*δs[k][l]/D[k][l];
				}
			}
			double h = (Fδ - G*Sδ)/(G*Dδ - Fδ*Sδ/Ss);
			
			for (int i = 0; i < n; i ++)
				for (int j = 0; j < m; j ++)
					g[i][j] = Math.max(0, g[i][j] + h/2*δg[i][j]);
			
			scorePrev = score;
			double L = 0;
			for (int k = 0; k < n; k ++)
				for (int l = 0; l < m; l ++)
					L += -1/2.*Math.pow(F[k][l] - G*s[k][l], 2)/D[k][l];
			double S = 0;
			for (int i = 0; i < n; i ++)
				for (int j = 0; j < m; j ++)
					if (g[i][j] > 1e-100)
						S += -g[i][j]*Math.log(g[i][j]);
			score = L + alpha*S;
			
			iter ++;
		} while (iter < 6 || score > scorePrev);
		
		double[][] s = gelfgatConvolve(P, g);
		double ΣFs = 0, Σss = 0;
		
		for (int k = 0; k < n; k ++) {
			for (int l = 0; l < m; l ++) {
				ΣFs += F[k][l]*s[k][l]/D[k][l]; // finalize the value of G
				Σss += s[k][l]*s[k][l]/D[k][l];
			}
		}
		G = ΣFs/Σss;
		
		double[][] source = new double[n][m];
		for (int i = 0; i < n; i ++)
			for (int j = 0; j < m; j ++)
				source[i][j] = G*g[i][j];
		return source; // finally, return the answer
	}


	private static double[][] gelfgatConvolve(double[][] P, double[][] g) {
		int n = g.length, m = g[0].length;
		double[][] s = new double[n][m];
		for (int i = 0; i < n; i ++)
			for (int j = 0; j < m; j ++)
				for (int k = 0; k < n; k ++)
					for (int l = 0; l < m; l ++)
						s[k][l] += P[m*k+l][m*i+j]*g[i][j];
		return s;
	}
	
	
	/**
	 * A two-dimensional array of numbers.
	 * 
	 * @author Justin Kunimune
	 */
	private static class Matrix {
		
		/**
		 * the array the contains the values of this Matrix.
		 */
		protected final double[][] values;
		
		/**
		 * Instantiate an nxm matrix with all zeroes.
		 * @param n - The height
		 * @param m - The width
		 */
		public Matrix(int n, int m) {
			this.values = new double[n][m];
		}

		/**
		 * Instantiate a matrix based on an existing 2D array
		 * @param values - The 2D array to load.
		 */
		public Matrix(double[][] values) {
			this.values = values;
		}
		
		/**
		 * Compute the transpose.
		 * @return the transpose of this
		 */
		public Matrix T() {
			Matrix tp = new Matrix(this.getM(), this.getN());
			for (int i = 0; i < tp.getN(); i ++)
				for (int j = 0; j < tp.getM(); j ++)
					tp.set(i, j, this.get(j, i));
			return tp;
		}
		
		/**
		 * Multiply this Matrix by an other matrix.
		 * @param that the other matrix
		 * @return the product
		 */
		public Matrix times(Matrix that) {
			if (this.getM() != that.getN())
				throw new IllegalArgumentException("Cannot multiply these matrices. Dimensions "+this.getN()+"x"+this.getM()+" and "+that.getN()+"x"+that.getM()+" do not agree.");
			Matrix prod = new Matrix(this.getN(), that.getM());
			for (int i = 0; i < this.getN(); i ++)
				for (int j = 0; j < that.getM(); j ++)
					for (int k = 0; k < this.getM(); k ++)
						prod.set(i, j, prod.get(i, j) + this.get(i, k)*that.get(k, j));
			return prod;
		}
		
		/**
		 * Multiply this Matrix by a scalar.
		 * @param a - The factor
		 * @return the product
		 */
		public Matrix times(double a) {
			Matrix product = new Matrix(this.getN(), this.getM());
			for (int i = 0; i < this.getN(); i ++)
				for (int j = 0; j < this.getM(); j ++)
					product.set(i, j, this.get(i, j) * a);
			return product;
		}
		
		/**
		 * A special kind of multiplication for column vectors.
		 * @param that - The Matrix to dot this with
		 * @return the dot product of this and that
		 */
		public double dot(Matrix that) {
			if (this.getM() != 1 || that.getM() != 1)
				throw new IllegalArgumentException("Cannot compute this dot product. Matrices "+this.getN()+"x"+this.getM()+" and "+that.getN()+"x"+that.getM()+" are not both column vectors.");
			if (this.getN() != that.getN())
				throw new IllegalArgumentException("Cannot compute this dot product. Vector lengths "+this.getN()+" and "+that.getN()+" do not match.");
			double product = 0;
			for (int i = 0; i < this.getN(); i ++)
				product += this.get(i, 0)*that.get(i, 0);
			return product;
		}
		
		/**
		 * Divide this Matrix by a scalar.
		 * @param a - The divisor
		 * @return the quotient
		 */
		public Matrix over(double a) {
			return this.times(1./a);
		}
		
		/**
		 * Add two Matrices.
		 * @param that - The addend
		 * @return the sum of this and that
		 */
		public Matrix plus(Matrix that) {
			if (this.getN() != that.getN() || this.getM() != that.getM())
				throw new IllegalArgumentException("Cannot add these matrices. Dimensions "+this.getN()+"x"+this.getM()+" and "+that.getN()+"x"+that.getM()+" do not match.");
			Matrix sum = new Matrix(this.getN(), this.getM());
			for (int i = 0; i < this.getN(); i ++)
				for (int j = 0; j < this.getM(); j ++)
					sum.set(i, j, this.get(i, j) + that.get(i, j));
			return sum;
		}
		
		/**
		 * Subtract two Matrices.
		 * @param that - The subtractee
		 * @return the difference of this and that
		 */
		public Matrix minus(Matrix that) {
			return this.plus(that.times(-1));
		}
		
		public Matrix inv() {
			if (this.getN() != this.getM())
				throw new IllegalArgumentException("Cannot invert a non-square matrix.");
			double[][] valuesClone = new double[this.getN()][this.getN()];
			for (int i = 0; i < this.getN(); i ++)
				for (int j = 0; j < this.getN(); j ++)
					valuesClone[i][j] = this.values[i][j];
			return new Matrix(invert(valuesClone));
		}
		
		/**
		 * L2 norm
		 * @return the sqrt of the sum of squares
		 */
		public double norm() {
			double s = 0;
			for (double[] row: values)
				for (double x: row)
					s += x*x;
			return Math.sqrt(s);
		}
		
		/**
		 * Compute the height.
		 * @return the height of this
		 */
		public int getN() {
			return this.values.length;
		}
		
		/**
		 * Compute the width.
		 * @return the width of this
		 */
		public int getM() {
			if (this.values.length > 0)
				return this.values[0].length;
			else
				return 0;
		}
		
		/**
		 * Extract a single scalar value.
		 * @return the value this_{i,j}
		 */
		public double get(int i, int j) {
			return this.values[i][j];
		}
		
		/**
		 * Extract a single row as a column vector.
		 * @return the value this_{i,j}
		 */
		public Matrix getRow(int i) {
			Matrix row = new Matrix(this.getM(), 1);
			for (int j = 0; j < this.getM(); j ++)
				row.set(j, 0, this.get(i, j));
			return row;
		}
		
		/**
		 * Set a single scalar value.
		 */
		public void set(int i, int j, double a) {
			this.values[i][j] = a;
		}
		
		public boolean equals(Object that) {
			if (!(that instanceof Matrix))
				return false;
			if (this.getN() != ((Matrix)that).getN() || this.getM() != ((Matrix)that).getM())
				return false;
			for (int i = 0; i < this.getN(); i ++)
				for (int j = 0; j < this.getM(); j ++)
					if (this.get(i, j) != ((Matrix)that).get(i, j))
						return false;
			return true;
		}
		
		public String toString() {
			StringBuilder str = new StringBuilder("Matrix(" + this.getN() + ", " + this.getM() + ", [\n");
			for (int i = 0; i < this.getN(); i ++) {
				str.append("  [");
				for (int j = 0; j < this.getM(); j ++)
					if (Math.abs(this.get(i, j)) >= 10000 || Math.abs(this.get(i, j)) < 0.1)
						str.append(String.format("%10.3e, ", this.get(i, j)));
					else
						str.append(String.format("% 10.4f, ", this.get(i, j)));
				str.append("],\n");
			}
			return str + "])";
		}
	}
	
	
	/**
	 * copied from https://www.sanfoundry.com/java-program-find-inverse-matrix/
	 */
	private static double[][] invert(double[][] a) {
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
	
	
//	/**
//	 * @param args
//	 */
//	public static void main(String[] args) {
////		Function<double[], Double> simionescu = (v) -> {
////			double x = v[0], y = v[1];
////			double z;
////			if (Math.hypot(x, y) > 1 + .2*Math.cos(8*Math.atan2(x, y)))
////				z = Double.POSITIVE_INFINITY;
////			else
////				z = .1*x*y;
////			System.out.printf("[%.4f, %.4f, %.4f],\n", x, y, z);
////			return z;
////		};
////		List<Function<double[], Double>> simionescuGradient = new ArrayList<Function<double[], Double>>(Arrays.asList(
////				(v) -> {
////					double y = v[1];
////					return .1*y;
////				},
////				(v) -> {
////					double x = v[0];
////					return .1*x;
////				}));
////		Function<double[], Double> himmelblau = (v) -> {
////			double x = v[0], y = v[1];
////			double z = Math.pow(x*x + y - 11, 2) + Math.pow(x + y*y - 7, 2) + 1;
////			System.out.printf("[%.4f, %.4f, %.4f],\n", x, y, z);
////			return z;
////		};
////		Function<double[], Double> ellipse = (v) -> {
////			double x = v[0], y = v[1];
////			double z = 2*Math.sqrt(1 + x*x - 1.5*x*y + y*y);
//////			System.out.printf("[%.4f, %.4f, %.4f],\n", x, y, z);
////			return z;
////		};
////		Function<double[], double[]> ellipseGrad = (v) -> {
////			double x = v[0], y = v[1];
////			double z = ellipse.apply(v);
////			return new double[] {
////				(2*x - 1.5*y)/(z/2),
////				(2*y - 1.5*x)/(z/2),
////			};
////		};
//		
////		System.out.println(Arrays.toString(minimizeNelderMead(
////				simionescu, new double[] {.5,.5}, 1e-8)));
////		System.out.println(Arrays.toString(minimizeCoordinateDescent(
////				simionescu, simionescuGradient, new double[] {.5, .5}, new double[] {1,1}, 1e-8)));
////		System.out.println(Arrays.toString(minimizeCoordinateDescent(
////				himmelblau, new double[] {0,0}, new double[] {1,1}, 1e-8)));
////		System.out.println(Arrays.toString(minimizeLBFGS(
////				himmelblau, new double[] {.5, .5}, 1e-8)));
////		System.out.println(Arrays.toString(minimizeLBFGSB(
////				ellipse, new double[] {1.5,0}, new double[] {.5,-1.5}, new double[] {2,1.5}, 1e-8)));
//	}

}
