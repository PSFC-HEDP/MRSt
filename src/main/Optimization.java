/**
 * 
 */
package main;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;
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
		if (grad >= 0)
			throw new IllegalArgumentException("Initial step must be downhill."); // if the gradient here is naught, there's absolutely noting we can do
		
		final double c = 0.4, τ = 0.5;
		
		double fx0 = func.apply(x0);
		double x = x0 + stepSize; // take an initial downhill step
		double fx = func.apply(x);
		while (fx - fx0 > (x - x0)*c*grad && Math.abs(x - x0) > x*1e-15) { // if the function didn't decrease enough (or we hit roundoff error)
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
	 * @param func0 the function at that point (I'm sure you already have it computed)
	 * @param grad0 the slope at that point (I'm not computing it myself)
	 * @param step0 the initial step length guess
	 * @param step1 the maximum allowed step length
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
		if (δf0 >= 0)
			throw new IllegalArgumentException("Initial step must be downhill.");// if the gradient here is naught, there's absolutely noting we can do
		
		System.out.println("wolf ");
		System.out.println("["+x0+", "+f0+", "+δf0+", "+x0+", "+(x0+stepMax)+"],");
		final double α = 1e-4, β = 0.9;
		
		double med = x0 + step0; // take an initial downhill step
		double min = x0, max = x0 + 2*stepMax;
		double truMax = x0 + stepMax;
		while (true) {
			if (max - min < 1e-15) // make sure it doesn't get stuck with an impossible function
				throw new RuntimeException("Could not find suitable minimum.");
			med = Math.min(med, truMax); // enforce that it not go past its true maximum
			double f = func.apply(med);
			System.out.println("["+med+", "+f+", "+grad.apply(med)+", "+min+", "+max+"],");
			if (f > f0 + α*(med - x0)*δf0) { // if the decrease condition is not met
				max = med; // we need to go closer
				med = (max + min)/2;
				continue;
			}
			if (max - min <= 1e-10*step0) { // if this has become quite tight
				return min; // there's most likely a discontinuity
			}
			double δf = grad.apply(med);
			if (Math.abs(δf) > β*Math.abs(δf0)) { // if the curvature condition is not met
				if (δf > 0) {
					max = med; // we might need to go closer
					med = (max + min)/2;
				}
				else if (med == truMax) {
					System.out.println("this boundary is likely the best we can do");
					return truMax; // this might be the best we can get
				}
				else {
					min = med; // but most likely the local min is farther out
					med = Math.min((max + min)/2, 2.7*min);
				}
				continue;
			}
			
			System.out.println("awøøøøø");
			return med; // if both are met, we're done here
		}
	}
	
	
	/**
	 * perform Nelder-Mead unconstrained optimization, making a SWAG for the initial simplex
	 * scale.
	 * @param func the function to minimize
	 * @param x0 the initial guess of the minimizing parameter vector
	 * @param scale an array of the same size as x0 giving the general length scale of its
	 * derivatives, to kick off the simplex.
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
			throw new IllegalArgumentException("Initial guess yielded bunk value");
		
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
					continue;
				}
				else {
					x[iWorst] = xR;
					fx[iWorst] = fxR;
					continue;
				}
			}
			else if (fxR < fx[iNext]) { // if this is better than the second worst
				x[iWorst] = xR;
				fx[iWorst] = fxR;
				continue;
			}
			else if (fxR < fx[iWorst]) { // if this has just a little bit of redeeming quality
				Matrix xS = xC.plus(xR.minus(xC).times(ρ));
				double fxS = f.apply(xS);
				
				if (fxS <= fxR) {
					x[iWorst] = xS;
					fx[iWorst] = fxS;
					continue;
				}
				else {
					for (int i = 0; i < x.length; i ++) {
						if (i != iBest) {
							x[i] = x[iBest].plus(x[i].minus(x[iBest]).times(σ)); // move all vertices inward
							fx[i] = f.apply(x[i]);
						}
					}
					continue;
				}
			}
			else {
				Matrix xS = xC.plus(x[iWorst].minus(xC).times(ρ)); // compute the contracted point
				double fxS = f.apply(xS);
				
				if (fxS < fx[iWorst]) {
					x[iWorst] = xS;
					fx[iWorst] = fxS;
					continue;
				}
				else {
					for (int i = 0; i < x.length; i ++) {
						if (i != iBest) {
							x[i] = x[iBest].plus(x[i].minus(x[iBest]).times(σ)); // move all vertices inward
							fx[i] = f.apply(x[i]);
						}
					}
					continue;
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
	 * @param func the function to minimize
	 * @param x0 the initial guess of the minimizing parameter vector
	 * @param lower the lowest allowable values of each of the elements of x
	 * @param upper the greatest allowable values of each of the elements of x
	 * @param tol the relative tolerance for error
	 * @return the vector x that minimizes f(x)
	 */
	public static double[] minimizeLBFGSB(
			Function<double[], Double> func, double[] x0, double[] lower, double[] upper, double tol) {
		final double ds = 1e-8;
		Function<double[], double[]> gradient = (x) -> { // finite difference gradient:
			double y0 = func.apply(x);
			double[] dydx = new double[x.length];
			for (int i = 0; i < x.length; i ++) { // compute each dimension individually
				double dx = (x[i] + ds >= upper[i]) ? -ds : ds; // chose a direction to avoid the bounds
				x[i] += dx; // perturb
				dydx[i] = (func.apply(x) - y0)/dx; // measure
				x[i] -= dx; // unperturb
			}
			return dydx;
		};
		BiFunction<double[], double[], Double> linGradient = (x, v) -> { // finite difference linear gradient:
			double vMag = 0; // this requires taking v's magnitude (to maintain the finite difference)
			for (int i = 0; i < v.length; i ++)
				vMag += v[i]*v[i];
			vMag = Math.sqrt(vMag);
			double y0 = func.apply(x);
			double dx = ds; // the step direction must have the same sign for all components now
			for (int i = 0; i < x.length; i ++)
				if (x[i] + v[i]/vMag*dx > upper[i] || x[i] + v[i]/vMag*dx < lower[i])
					dx *= -1; // but we can still flip it if one of the components is near a bound
			for (int i = 0; i < x.length; i ++)
				x[i] += v[i]/vMag*dx; // do all the dimensionless perturbations
			double dydx = (func.apply(x) - y0)/dx; // measure
			for (int i = 0; i < x.length; i ++)
				x[i] -= v[i]/vMag*dx; // then put them all back
			return dydx*vMag; // remember to scale by vMag since this is a dot product
		};
		return minimizeLBFGSB(func, gradient, linGradient, x0, lower, upper, tol);
	}
	
	
	/**
	 * perform L-BFGS-B constrained optimization with a provided gradient.
	 * <p>
	 * Richard H. Byrd, Peihuang Lu, Jorge Nocedal, & Ciyou Zhu (1994). "A limited memory
	 *   algorithm for bound constrained optimization." <i>Northwestern University Department
	 *   of Electrical Engineering and Computer Science.</i> Technical Report NAM-08.
	 * @param func the function to minimize
	 * @param the gradient of that function
	 * @param x0 the initial guess of the minimizing parameter vector
	 * @param lower the lowest allowable values of each of the elements of x
	 * @param upper the greatest allowable values of each of the elements of x
	 * @param tol the relative tolerance for error
	 * @return the vector x that minimizes f(x)
	 */
	public static double[] minimizeLBFGSB(
			Function<double[], Double> func, Function<double[], double[]> grad, double[] x0,
			double[] lower, double[] upper, double tol) {
		return minimizeLBFGSB(
				func, grad, (x, v) -> {
					double[] g = grad.apply(x);
					double gv = 0;
					for (int i = 0; i < v.length; i ++)
						gv += g[i]*v[i];
					return gv;
				}, x0, lower, upper, tol);
	}
	
	/**
	 * Perform L-BFGS-B constrained optimization with a provided multidimensional gradient. In
	 * addition, since it may be much faster to compute the gradient in a single direction,
	 * this function also accepts a separate linear gradient function.
	 * @param func the function to minimize
	 * @param grad the gradient of that function
	 * @param linGrad the gradient of that function dotted with the second input
	 * @param x0 the initial guess of the minimizing parameter vector
	 * @param lower the lowest allowable values of each of the elements of x
	 * @param upper the greatest allowable values of each of the elements of x
	 * @param tol the relative tolerance for error
	 * @return the vector x that minimizes f(x)
	 */
	public static double[] minimizeLBFGSB(
			Function<double[], Double> func, Function<double[], double[]> grad, BiFunction<double[], double[], Double> linGrad,
			double[] x0, double[] lower, double[] upper, double tol) {
		if (lower.length != x0.length || upper.length != x0.length)
			throw new IllegalArgumentException("Bound lengths don't match");
		for (int i = 0; i < x0.length; i ++)
			if (lower[i] > x0[i] || x0[i] > upper[i])
				throw new IllegalArgumentException("Upper bounds must be greater than lower bounds");
		
		final int mMax = 6;
		final int n = x0.length;
		
		Function<Matrix, Double> funcMat = (mat) -> func.apply(mat.T().values[0]); // change these from double[] to Matrix things
		Function<Matrix, Matrix> gradMat = (mat) -> new Matrix(new double[][] {grad.apply(mat.T().values[0])}).T();
		BiFunction<Matrix, Matrix, Double> linGradMat = (mat, v) -> linGrad.apply(mat.T().values[0], v.T().values[0]);
		
		LinkedList<Matrix> yHist = new LinkedList<Matrix>();
		LinkedList<Matrix> sHist = new LinkedList<Matrix>();
		double θ = 1;
		Matrix xk = new Matrix(new double[][] {x0.clone()}).T(); // the current best guess (column matrix)
		Matrix gk = gradMat.apply(xk); // the current value of g from the L-BFGS algorithm
		double fxk = funcMat.apply(xk);
		if (!Double.isFinite(fxk))
			throw new IllegalArgumentException("Initial guess yielded bunk value");
		
		while (true) {
			System.out.println(Arrays.toString(xk.T().values[0]));
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
					for (int j = 0; j < m; j ++)
						if (i > j)
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
				Matrix B0 = new Matrix(n, n);
				for (int i = 0; i < n; i ++)
					B0.set(i, i, θ);;
				Matrix Bk = B0.minus(Wk.times(Mk.times(Wk.T())));
				
				for (int i = 0; i < 10000; i ++) {
					Matrix u = new Matrix(n, 1);
					for (int j = 0; j < n; j ++)
						u.set(j, 0, 2*Math.random() - 1);
					assert u.dot(Bk.times(u)) > 0: Bk;
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
				LinkedList<Integer> breakpointOrder = new LinkedList<Integer>();
				for (int i = 0; i < n; i ++) // in terms of the original paper, this list contains the order in which indices are removed from F
					if (breakpoints[i] > 0)
						breakpointOrder.add(i);
				breakpointOrder.sort((iA, iB) -> (int)Math.signum(breakpoints[iA] - breakpoints[iB]));
//				System.out.println(breakpointOrder);
				
				Matrix p = Wk.T().times(d); // first look for minimum in first segment
				Matrix c = new Matrix(2*m, 1);
				double dfdt = gk.dot(d);
				double d2fdt2 = -θ*dfdt - p.dot(Mk.times(p));
				double Δtmin = -dfdt/d2fdt2;
				double told = 0;
				int b = breakpointOrder.pop();
				double t = breakpoints[b];
				double Δt = t;
//				System.out.println("f′ = "+dfdt+" = "+d.dot(gk));
//				System.out.println("f″ = "+d2fdt2+" = "+d.dot(Bk.times(d)));
//				System.out.println("["+told+", "+t+", "+(gk.dot(xC.minus(xk))+0.5*xC.minus(xk).dot(Bk.times(xC.minus(xk))))+", "+dfdt+", "+d2fdt2+"],");
				while (Δtmin >= Δt) { // then check all subsequent segments
//					System.out.println("your minimum is in another interval!");
					double xCb = (d.get(b, 0) > 0) ? upper[b] : lower[b];
					double zb = xCb - xk.get(b, 0);
					double gb = gk.get(b, 0);
					Matrix wb = Wk.getRow(b);
					c = c.plus(p.times(Δt));
					dfdt = dfdt + Δt*d2fdt2 + gb*gb + θ*gb*zb - gb*wb.dot(Mk.times(c));
					d2fdt2 = d2fdt2 - θ*gb*gb - 2*gb*wb.dot(Mk.times(p)) - gb*gb*wb.dot(Mk.times(wb));
					assert d2fdt2 > 0 : d2fdt2;
					p = p.plus(wb.times(gb));
					d.set(b, 0, 0);
					Δtmin = -dfdt/d2fdt2;
					told = t;
					b = breakpointOrder.pop();
					t = breakpoints[b];
					Δt = t - told;
//					System.out.println("f′ = "+dfdt+" = "+d.dot(gk.plus(Bk.times(xC.minus(xk)))));
//					System.out.println("f″ = "+d2fdt2+" = "+d.dot(Bk.times(d)));
//					System.out.println("["+told+", "+t+", "+(gk.dot(xC.minus(xk))+0.5*xC.minus(xk).dot(Bk.times(xC.minus(xk))))+", "+dfdt+", "+d2fdt2+"],");
				}
				Δtmin = Math.max(Δtmin, 0);
				told = told + Δtmin;
//				System.out.println("Looks like the min ended up being at "+told);
				c = c.plus(p.times(Δtmin));
				List<Integer> F = new ArrayList<Integer>(breakpointOrder.size()+1);
				Matrix xC = new Matrix(n, 1);
				for (int i = 0; i < n; i ++) { // I'm setting xC here not how it's done in the paper, because the paper version is _totally_ wrong for this part
					if (breakpoints[i] <= told)
						xC.set(i, 0, (gk.get(i, 0) < 0) ? upper[i] : lower[i]);
					else
						xC.set(i, 0, xk.get(i, 0) - told*gk.get(i, 0));
					if (xC.get(i, 0) != lower[i] && xC.get(i, 0) != upper[i])
						F.add(i);
				}
//				System.out.println(xk.T());
//				System.out.println(xC.T());
//				System.out.println(gk.dot(xC.minus(xk)) + 0.5*xC.minus(xk).dot(Bk.times(xC.minus(xk))));
//				System.out.println(c.T());
//				System.out.println(Wk.T().times(xC.minus(xk)).T());
				System.out.println(F.size()+"/"+n);
				assert xC.equals(P(xC, lower, upper));
				assert xC.minus(xk).dot(gk) < 0;
				assert xC.minus(xk).dot(Bk.times(xC.minus(xk))) > 0;
				assert fxk + gk.dot(xC.minus(xk)) + 0.5*xC.minus(xk).dot(Bk.times(xC.minus(xk))) < fxk;
				
				if (F.size() >= 1) {
					System.out.println(θ);
					Matrix Zk = new Matrix(n, F.size()); // STEP 3: find an approximate bound minimum (direct primal method)
					for (int f = 0; f < F.size(); f ++)
						Zk.set(F.get(f), f, 1);
					assert Zk.T().times(Zk).inv().equals(Zk.T().times(Zk));
					Matrix BHatk = Zk.T().times(Bk.times(Zk));
					Matrix rHatC = Zk.T().times(gk.plus(xC.minus(xk).times(θ)).minus(Wk.times(Mk.times(c))));
//					System.out.println(Zk.T().times(gk.plus(Bk.times(xC.minus(xk)))).T());
//					System.out.println(rHatC.T());
					Matrix v = Wk.T().times(Zk.times(rHatC));
					v = Mk.times(v);
					Matrix N = Mk.times(Wk.T().times(Zk.times(Zk.T().times(Wk)))).over(-θ);
					assert N.getN() == 2*m && N.getM() == 2*m;
					for (int i = 0; i < 2*m; i ++)
						N.set(i, i, 1 + N.get(i, i));
					v = N.inv().times(v);
					Matrix dHatU = rHatC.over(θ).plus(Zk.T().times(Wk.times(v)).over(θ*θ)).times(-1); // the paper has a sign error
					for (int i = 0; i < F.size(); i ++)
						System.out.print(dHatU.get(i, 0) / BHatk.inv().times(rHatC).times(-1).get(i, 0) + " ");
					System.out.println();
					Matrix dU = Zk.times(dHatU);
					double αStar = 1;
					for (int i = 0; i < n; i ++) {
						if (xC.get(i, 0) + dU.get(i, 0) > upper[i]) // the paper has an error here; you need to multiply dU by Zk before comparing to u and l
							αStar = Math.min(αStar, (upper[i] - xC.get(i, 0))/dU.get(i, 0));
						else if (xC.get(i, 0) + dU.get(i, 0) < lower[i])
							αStar = Math.min(αStar, (lower[i] - xC.get(i, 0)/dU.get(i, 0)));
					}
					System.out.println(αStar);
					assert αStar > 0 && αStar <= 1 : αStar;
					Matrix xBar = xC.plus(dU.times(αStar));
					xBar = P(xBar, lower, upper); // strictly speaking I shouldn't need this, but roundoff
	//				System.out.println(gk.T());
					System.out.println("[0, "+(1/αStar)+", "+(fxk + gk.dot(xC.minus(xk)) + xC.minus(xk).dot(Bk.times(xC.minus(xk))))+", "+(dU.times(αStar).dot(gk.plus(Bk.times(xC.minus(xk)))))+", "+(dU.times(αStar).dot(Bk.times(dU.times(αStar))))+"]");
					if (xBar.minus(xk).dot(Bk.times(xBar.minus(xk))) <= 0) {
						System.out.println("not positive definite!!");
						System.out.println(xBar.minus(xk).T());
						System.out.println(Arrays.deepToString(Bk.values));
						System.out.println(xBar.minus(xk).dot(Bk.times(xBar.minus(xk))));
						Bk = B0;
						for (int j = 0; j < m; j ++) {
							Matrix s = sHist.get(j);
							Matrix y = yHist.get(j);
							Bk = Bk.minus(Bk.times(s.times(s.T().times(Bk))).over(s.dot(Bk.times(s)))).plus(y.times(y.T()).over(y.dot(s)));
						}
						System.out.println(Arrays.deepToString(Bk.values));
						System.out.println(xBar.minus(xk).dot(Bk.times(xBar.minus(xk))));
						assert false;
					}
					
					if (gk.dot(xBar.minus(xk)) + 0.5*xBar.minus(xk).dot(Bk.times(xBar.minus(xk))) > gk.dot(xC.minus(xk)) + 0.5*xC.minus(xk).dot(Bk.times(xC.minus(xk)))) {
//						System.out.println(Arrays.deepToString(Bk.values));
//						Bk = B0;
//						for (int j = 0; j < m; j ++) {
//							Matrix s = sHist.get(j);
//							Matrix y = yHist.get(j);
//							Bk = Bk.minus(Bk.times(s.times(s.T().times(Bk))).over(s.dot(Bk.times(s)))).plus(y.times(y.T()).over(y.dot(s)));
//						}
//						System.out.println(Arrays.deepToString(Bk.values));
//						System.out.println(Arrays.deepToString(BHatk.values));
						System.out.println(Arrays.deepToString(Zk.T().times(Bk.times(Zk)).values));
						Matrix BHatL = Zk.T().times(Wk).times(Mk.times(Wk.T().times(Zk))).times(-1);
						for (int i = 0; i < BHatL.getN(); i ++)
							BHatL.set(i, i, θ + BHatL.get(i, i));
						System.out.println(Arrays.deepToString(BHatL.values));
						System.out.println("invert!");
						System.out.println(Arrays.deepToString(BHatk.inv().values));
						Matrix BHatinv = new Matrix(F.size(), F.size());
						for (int i = 0; i < F.size(); i ++)
							BHatinv.set(i, i, 1/θ);
						Matrix inner = new Matrix(2*m, 2*m);
						for (int i = 0; i < 2*m; i ++)
							inner.set(i, i, 1);
						inner = inner.minus(Mk.times(Wk.T().times(Zk.times(Zk.T().times(Wk)))).over(θ));
						BHatinv = BHatinv.plus(Zk.T().times(Wk.times(inner.inv().times(Mk.times(Wk.T().times(Zk))))).over(θ*θ));
						System.out.println(Arrays.deepToString(BHatinv.values));
						System.out.println("just to check...");
						System.out.println(Arrays.deepToString(inner.times(inner.inv()).values));
						System.out.println(Arrays.deepToString(inner.inv().times(inner).values));
						System.out.println("ANd then some other stuff me lew jana esa galti");
						System.out.println(BHatk.inv().times(rHatC).times(-1).T());
						System.out.println(dHatU.T());
						System.out.println(" is the step and the new gradient should be ");
						System.out.println(Zk.times(rHatC.plus(BHatk.times(BHatk.inv().times(rHatC).times(-1)))).T());
						System.out.println(Zk.times(rHatC.plus(BHatk.times(dHatU))).T());
						System.out.println(Zk.times(rHatC.plus(Zk.T().times(Bk.times(Zk.times(dHatU))))).T());
						System.out.println(Zk.times(rHatC).plus(Zk.times(Zk.T().times(Bk.times(Zk.times(dHatU))))).T());
						System.out.println(Zk.times(Zk.T().times(gk.plus(Bk.times(xC.minus(xk))))).plus(Zk.times(Zk.T().times(Bk.times(Zk.times(dHatU))))).T());
						System.out.println(gk.plus(Bk.times(xC.minus(xk))).plus(Bk.times(xBar.minus(xC))).T());
						System.out.println(gk.plus(Bk.times(xC.minus(xk))).plus(Bk.times(xBar.minus(xC))).T());
						System.out.println(gk.plus(Bk.times(xBar.minus(xk))).T());
						System.out.println("in conclusion...");
						Matrix A = new Matrix(F.size(), F.size());
						for (int i = 0; i < F.size(); i ++)
							A.set(i, i, θ);
						Matrix U = Zk.T().times(Wk);
						Matrix C = Mk.times(-1);
						Matrix V = Wk.T().times(Zk);
						System.out.println("A = "+A);
						System.out.println("U = "+U);
						System.out.println("C = "+C);
						System.out.println("V = "+V);
						System.out.println("CV = "+C.times(V));
						double sum = 0;
						for (int i = 0; i < V.getN(); i ++) {
							sum += C.get(0, i)*V.get(i, 0);
							System.out.print(C.get(0, i)+"*"+V.get(i, 0)+" + ");
						}
						System.out.println("= "+sum);
						System.out.println("UCV = "+U.times(C.times(V)));
						System.out.println("A + UCV = "+A.plus(U.times(C.times(V))));
						System.out.println("(A + UCV)^-1 = "+A.plus(U.times(C.times(V))).inv());
						System.out.println("A^-1 - A^-1U(C^-1 + VA^-1U)^-1VA^-1 = "+BHatinv);
					}
	//				System.out.println("That looks like a minimum.");
					//				for (int i = 0; i < n; i ++) {
	//					assert xBar.get(i, 0) >= lower[i] : xBar.get(i, 0)+"<"+lower[i]; // if this is violated it's only because of roundoff
	//					assert xBar.get(i, 0) <= upper[i] : xBar.get(i, 0)+">"+upper[i];
	//				}
					for (int i = 0; i < n; i ++)
						assert F.contains(i) || xBar.get(i, 0) == xC.get(i, 0);
					assert gk.dot(xBar.minus(xk)) + 0.5*xBar.minus(xk).dot(Bk.times(xBar.minus(xk))) <= gk.dot(xC.minus(xk)) + 0.5*xC.minus(xk).dot(Bk.times(xC.minus(xk))) :
						(fxk+gk.dot(xBar.minus(xk)) + 0.5*xBar.minus(xk).dot(Bk.times(xBar.minus(xk))))+" > "+(fxk+gk.dot(xC.minus(xk)) + 0.5*xC.minus(xk).dot(Bk.times(xC.minus(xk))));
					dk = xBar.minus(xk);
				}
				else {
					dk = xC.minus(xk);
				}
			}
			else {
				dk = P(xk.minus(gk), lower, upper).minus(xk);
			}
			assert gk.dot(dk) <= 0;
			
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
			final Matrix Xk = xk;
			double λk = minimizeWolfe(
					(λ) -> {
						return funcMat.apply(P(Xk.plus(dk.times(λ)), lower, upper)); // strictly speaking I shouldn't need to call P here, but there's roundoff stuff
					}, (λ) -> {
						return linGradMat.apply(P(Xk.plus(dk.times(λ)), lower, upper), dk);
					}, 0, fxk, gk.dot(dk), 1, λMax);
			xk = P(xk.plus(dk.times(λk)), lower, upper);
			
			Matrix gkp1 = gradMat.apply(xk);
			double fxkp1 = funcMat.apply(xk);
			assert fxkp1 <= fxk;
			if (sHist.size() > 1 && Math.abs((fxkp1 - fxk)/fxk) < tol) { // STEP 5: stop condition
				return xk.T().values[0]; // if we're into it and the energy isn't really changing, then we're done
			}
			
			Matrix sk = dk.times(λk); // STEP 6: save historical vector information
			Matrix yk = gkp1.minus(gk);
			if (yk.dot(sk) > 1e-15*yk.norm()) {
				yHist.addLast(yk);
				sHist.addLast(sk);
				θ = yk.dot(yk)/yk.dot(sk);
			}
			else
				System.out.println("SKIPPING");
			if (yHist.size() > mMax) {
				yHist.removeFirst();
				sHist.removeFirst();
			}
			gk = gkp1;
			fxk = fxkp1;
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
		 * Multiply this Matrix by a scalar.
		 * @param a - The factor
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
		 * @param i
		 * @param j
		 * @return the value this_{i,j}
		 */
		public double get(int i, int j) {
			return this.values[i][j];
		}
		
		/**
		 * Extract a single row as a column vector.
		 * @param i
		 * @param j
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
		 * @param i
		 * @param j
		 * @param a
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
			String str = "Matrix("+this.getN()+", "+this.getM()+", [\n";
			for (int i = 0; i < this.getN(); i ++) {
				str += "  [";
				for (int j = 0; j < this.getM(); j ++)
					if (Math.abs(this.get(i, j)) >= 10000 || Math.abs(this.get(i, j)) < 0.1)
						str += String.format("%10.3e, ", this.get(i, j));
					else
						str += String.format("% 10.4f, ", this.get(i, j));
				str += "],\n";
			}
			return str + "])";
		}
	}
	
	
	/**
	 * copied from https://www.sanfoundry.com/java-program-find-inverse-matrix/
	 * @param a
	 * @return
	 */
	public static double[][] invert(double a[][]) {
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
	 * @param args
	 */
	public static void main(String[] args) {
		/*Function<double[], Double> simionescu = (v) -> {
			double x = v[0], y = v[1];
			double z;
			if (Math.hypot(x, y) > 1 + .2*Math.cos(8*Math.atan2(x, y)))
				z = Double.POSITIVE_INFINITY;
			else
				z = .1*x*y;
			System.out.printf("[%.4f, %.4f, %.4f],\n", x, y, z);
			return z;
		};
		List<Function<double[], Double>> simionescuGradient = new ArrayList<Function<double[], Double>>(Arrays.asList(
				(v) -> {
					double y = v[1];
					return .1*y;
				},
				(v) -> {
					double x = v[0];
					return .1*x;
				}));
//		Function<double[], Double> himmelblau = (v) -> {
//			double x = v[0], y = v[1];
//			double z = Math.pow(x*x + y - 11, 2) + Math.pow(x + y*y - 7, 2) + 1;
//			System.out.printf("[%.4f, %.4f, %.4f],\n", x, y, z);
//			return z;
//		};
		
//		System.out.println(Arrays.toString(minimizeNelderMead(
//				simionescu, new double[] {.5,.5}, 1e-8)));
		System.out.println(Arrays.toString(minimizeCoordinateDescent(
				simionescu, simionescuGradient, new double[] {.5, .5}, new double[] {1,1}, 1e-8)));
//		System.out.println(Arrays.toString(minimizeCoordinateDescent(
//				himmelblau, new double[] {0,0}, new double[] {1,1}, 1e-8)));
//		System.out.println(Arrays.toString(minimizeLBFGS(
//				himmelblau, new double[] {.5, .5}, 1e-8)));*/
		Matrix C = new Matrix(new double[][] {
  { 0.0206,  0.0025,  0.0018,  0.0061, -0.0010,  0.0000,  0.0929, -0.0201, -0.1133,  0.0319, -0.0055,  0.0026, },
  { 0.0025,  0.0006,  0.0017,  0.0099,  0.0002,  0.0000, -0.0011,  0.0101, -0.0220, -0.0325,  0.0070, -0.0005, },
  { 0.0018,  0.0017,  0.0134,  0.0503, -0.0043,  0.0000, -0.0026, -0.0014,  0.2858, -0.2738, -0.0095,  0.0110, },
  { 0.0061,  0.0099,  0.0503,  0.2675,  0.0007,  0.0000, -0.0088, -0.0087, -0.0164,  0.0843, -0.0054, -0.0018, },
  {-0.0010,  0.0002, -0.0043,  0.0007,  0.0095,  0.0000,  0.0012, -0.0004,  0.0033, -0.0034, -0.0003,  0.0002, },
  { 0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0020,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, },
  { 0.0929, -0.0011, -0.0026, -0.0088,  0.0012, -0.0000, -0.1076,  0.0215,  0.1358, -0.0315,  0.0051, -0.0030, },
  {-0.0201,  0.0101, -0.0014, -0.0087, -0.0004, -0.0000,  0.0215, -0.0146, -0.0025,  0.0399, -0.0083,  0.0011, },
  {-0.1133, -0.0220,  0.2858, -0.0164,  0.0033, -0.0000,  0.1358, -0.0025, -0.5084,  0.2394,  0.0191, -0.0084, },
  { 0.0319, -0.0325, -0.2738,  0.0843, -0.0034, -0.0000, -0.0315,  0.0399,  0.2394, -0.3992,  0.0167,  0.0086, },
  {-0.0055,  0.0070, -0.0095, -0.0054, -0.0003,  0.0000,  0.0051, -0.0083,  0.0191,  0.0167, -0.0056,  0.0008, },
  { 0.0026, -0.0005,  0.0110, -0.0018,  0.0002,  0.0000, -0.0030,  0.0011, -0.0084,  0.0086,  0.0008, -0.0005, },
});
		Matrix V = new Matrix(new double[][] {
  { 0.0087,  0.4235,  2.0722,  3.7689,  2.4753,  0.8018,  0.2256,  0.3623,  0.0146,  0.0102,  0.0218,  0.0029, -0.1979, -0.1746, -0.1513,  0.0029, -0.0931, -0.5937, -0.2692, -58379.1698,  4.5140,  3.8301,  0.4322,  0.0000, -0.0015, -0.0058,  0.0000,  0.0015,  0.0029,  0.0160,  0.0378,  9.3991, -0.0553, -0.1251, -0.0116, -0.0204, -1.0012, -0.9706, -0.6039,  0.2488, -1.9412, -1.8510,  3.9727, -93.6969,  10.4061,  82.4220,  9.6552, },
  { 0.0204,  0.9852,  5.0059,  9.2434,  5.2896,  0.6257, -0.1994,  0.6708, -0.0524,  0.0247,  0.0509,  0.0058, -0.5544, -0.4511, -0.4075, -0.2357, -0.4482,  1.3853,  7.8187, -178459.1179, -33.6673,  6.3956,  0.3332,  0.0000,  0.0000,  0.0029,  0.0015,  0.0233,  0.0655,  0.1513,  0.3769,  32.3183,  0.2648, -0.1863, -0.0204, -0.0029, -2.7270, -2.5742, -2.0678, -1.8626, -3.1316,  14.9244, -19.5781, -355.2967, -133.5779,  198.8272,  22.9935, },
  { 0.0058,  0.2023,  1.0419,  1.8903,  1.1758,  0.3944,  0.0713, -0.2328, -0.0160, -0.0044,  0.0044,  0.0000, -0.1120, -0.0917, -0.0946,  0.0320,  0.1572, -0.2663, -0.0364,  48210.6494,  7.0082,  2.3065,  0.2357,  0.0000,  0.0015,  0.0029,  0.0044,  0.0044,  0.0044, -0.0058, -0.0276, -9.2215, -0.0553, -0.0655, -0.0102, -0.0189, -0.5472, -0.5472, -0.3653,  0.3405, -0.5035, -1.6196,  2.2686,  97.2475,  21.5805,  46.6607,  5.1310, },
  { 0.0000,  0.2357,  1.2369,  2.2832,  1.3912,  0.4191,  0.0728,  0.0538, -0.0029,  0.0073,  0.0131,  0.0015, -0.1135, -0.1164, -0.1251, -0.0087,  0.1106, -0.1455,  0.5341, -8546.4257,  2.1522,  2.1260,  0.1717,  0.0000, -0.0029, -0.0058,  0.0000,  0.0029,  0.0087,  0.0131,  0.0189,  1.6516, -0.0029, -0.0597, -0.0058, -0.0175, -0.6607, -0.6650, -0.5108,  0.0742, -0.7276, -0.1179,  0.6548, -12.3400, -0.3143,  48.7387,  5.2838, },
  {-0.0073,  0.5661,  3.0341,  5.5705,  4.4398,  3.0777,  1.3781,  0.8178,  0.1251,  0.0306,  0.0044,  0.0058, -0.3245, -0.2750, -0.3783,  0.3027,  1.0783, -4.4180, -10.8164, -46822.1406,  65.9493,  8.0108,  1.1220,  0.0000,  0.0029,  0.0073, -0.0015, -0.0146, -0.0378, -0.1412, -0.3391,  9.4704, -0.5792, -0.2925, -0.0218, -0.1412, -1.3097, -1.5600, -0.8717,  3.5507, -3.8097, -27.7374,  45.2215,  31.6519,  204.8968,  104.4406,  10.2111, },
  { 0.0073,  1.2165,  6.9238,  13.4562,  8.2873,  2.4243,  0.7538,  6.2690,  0.1950,  0.1848,  0.0888,  0.0087, -0.8557, -0.7174, -0.9619, -0.7058,  0.2052,  0.4482,  7.9264, -1293538.0990, -53.7533,  3.2975, -0.5806,  0.0000, -0.0044, -0.0029,  0.0073,  0.0538,  0.1339,  0.2707,  0.8906,  292.9519,  0.6912, -0.1208, -0.0160, -0.1004, -3.9698, -4.2128, -4.5664, -4.0309, -6.6138,  18.1346, -13.6119, -2222.1750, -367.5785,  161.6660,  13.7865, },
  { 1043563714.0429,  802834531.1439,  604099010.8408,  592237731.4698,  583451944.3369,  485682168.0901,  312680725.1851,  141798379.5559,  38020440.7679,  1617192.4037, -1737535.1616, -9092.5871, -5378602.4234, -1620997.4555, -5802640.0997, -10865701.6850, -17432255.9541, -15206314.0272, -22612200.2138, -156986.4858,  79838831.2910,  577175756.8235,  145675469.5388,  0.0000, -81124.8025, -144895.5049, -268655.6016, -1040369.1299, -5046776.6778, -22136982.9365, -76043985.7419, -150674672.0235, -130756604.4074, -39408345.6798, -2943868.2247,  105090.5683,  2923323.4868,  2328.5089, -10927510.0290, -38798651.6428, -67067731.4772,  5075447.0747,  29937414.7868,  47574813.4415, -37692012.8998,  247538301.7130,  205745011.5903, },
  { 2747445836.3548,  2112490795.2248,  1572975390.3015,  1507467561.4374,  1478811197.4114,  1250451785.3174,  810399969.2844,  370016915.3320,  101410146.5153,  5805358.2260, -3857029.8171, -4752.6669, -11807020.4145, -4789605.9530, -12597100.5391, -18184533.3534, -29988424.2693, -50181132.5113,  6658558.3320,  918818.0153, -227930418.7722,  1510297728.2248,  380207579.9921,  0.0000, -193915.0601, -375119.8082, -689751.5674, -2635166.7690, -12876233.4872, -57060943.8115, -198692928.1440, -395822341.9247, -340011362.7276, -102638794.7294, -7754649.0320,  89761.3270, -1176843.8686, -3427244.7605, -18195601.4234, -79616427.6951, -211760782.5005, -492476681.5306, -1185332128.8983, -641725809.8664, -440750669.7632,  746452353.1259,  675585361.5751, },
  { 585802570.7046,  450801081.5457,  340951701.6443,  337760942.9375,  333203529.6142,  275313402.9558,  176683267.4856,  79849157.4758,  21198379.3793,  745275.0383, -1048078.2094, -7014.7069, -3268359.4095, -865445.3321, -3560028.6067, -7282052.5092, -11823983.5805, -8772290.6093, -18359456.8929, -127503.8874,  77665858.3258,  324920406.8085,  82094598.1843,  0.0000, -47551.8363, -81940.6025, -152054.7921, -591671.4079, -2863079.8607, -12504774.4343, -42747676.1733, -84679690.8072, -73821420.9037, -22226350.4427, -1651644.2320,  78898.7234,  2504162.7336,  295989.9238, -7568570.3479, -26132750.0685, -40579549.1316,  35047155.2953,  98391263.7015,  77831697.7316, -3382844.2294,  129373702.2385,  101778355.8554, },
  { 699764880.5713,  538437245.9210,  406202494.7781,  400088820.0568,  394039538.9250,  326750012.5567,  210023048.5670,  95055534.5594,  25395633.9078,  986610.8104, -1203967.0208, -7052.6213, -3760125.8333, -1079391.6674, -4112813.0096, -8161749.2437, -13601179.3061, -12763288.8592, -16179049.9443, -44297.7328,  48789591.0364,  387491674.8881,  97822182.0773,  0.0000, -55511.3143, -97405.1819, -179783.9342, -696188.4421, -3377933.0978, -14798103.5048, -50844926.4883, -101109107.0843, -87888754.8981, -26465196.5804, -1972640.4741,  82963.4068,  2351745.5586,  56197.4090, -8794050.9369, -32336773.9238, -59077419.4250, -16189144.4543, -26131351.0371,  7874745.2686, -48267995.6545,  160907667.9435,  130354596.3394, },
  { 1588708027.9021,  1224848449.5008,  956642034.0857,  1008456039.9308,  1002980441.5719,  793673735.2617,  499665566.6434,  221087443.7283,  54996013.4147, -801591.8684, -4110168.8455, -52399.1577, -13186366.2034, -1555202.2534, -14873083.5298, -40125172.0112, -66788885.2709, -25487106.5285, -150574061.6719, -2029933.2937,  805453311.5518,  897326087.1394,  228251835.1670,  0.0000, -164011.4371, -232825.7870, -435110.7699, -1743697.1434, -8311284.3067, -35363981.8607, -117165378.2539, -231353582.9221, -207568300.0875, -62116969.6945, -4464658.1308,  559189.4297,  21909288.9459,  6042484.8950, -44839546.5872, -142803085.2147, -148709442.2127,  692930547.1947,  1778384448.3467,  1153389127.9952,  332981390.5919,  183775448.9801,  37167726.2619, },
  { 4161995341.2438,  3202957479.2736,  2419190001.6248,  2383789089.5570,  2341240148.5607,  1937679436.3457,  1244577029.6009,  562419073.4729,  150584159.5437,  5606416.6386, -7186223.2169, -41396.4469, -22832038.8243, -6707355.9007, -25691330.7519, -53519343.0485, -97006579.3404, -121674664.9631, -60543901.8734,  575252.3821, -115866023.0827,  2304454393.8727,  581330110.7685,  0.0000, -331559.3896, -578085.9578, -1050458.8869, -4039889.0999, -19684478.4384, -86399918.4925, -299081753.2790, -601795727.3830, -522589088.5350, -157142663.2801, -11715435.1861,  535246.8992,  12984909.3896, -1485819.5987, -65907998.3168, -268615477.2574, -587109758.6340, -773892028.4818, -1808086598.3890, -896607227.4806, -889282699.6676,  954256501.9162,  765289578.9861, },
});
		System.out.println(C);
		System.out.println(V);
		System.out.println(C.times(V));
	}


}
