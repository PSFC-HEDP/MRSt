/**
 * 
 */
package main;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
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
	 * @param x0 the starting guess of x0
	 * @param stepSize the upper bound on the distance to the minimum (should be positive)
	 * @return x such that f(x) is better than f(x0)
	 */
	public static double minimizeBacktrack(
			Function<Double, Double> func, double grad, double x0, double stepSize) {
		if (grad == 0)  return x0; // if the gradient here is naught, there's absolutely noting we can do
		
		final double c = 0.4, τ = 0.5;
		
		double fx0 = func.apply(x0);
		if (grad > 0)  stepSize *= -1; // orient ourselves
		double x = x0 + stepSize; // take an initial downhill step
		double fx = func.apply(x);
		while (fx - fx0 > (x - x0)*c*grad && Math.abs(x - x0) > x*1e-15) { // if the function didn't decrease enough (or we hit roundoff error)
			x = x0 + τ*(x - x0); // backstep and try again
			fx = func.apply(x);
		}
		
		return x;
	}
	
	
	/**
	 * perform Nelder-Mead unconstrained optimization.
	 * @param func the function to minimize
	 * @param x0 the initial guess of the minimizing parameter vector
	 * @param tol the relative tolerance for error
	 * @return the x vector that minimizes f(x)
	 */
	public static double[] minimizeNelderMead(
			Function<double[], Double> func, double[] x0, double tol) {
		final double α = 1, γ = 2, ρ = 1/2., σ = 1/2.; // declare values
		
		Function<Matrix, Double> f = (mat) -> func.apply(mat.values[0]); // change these from double[] to Matrix things
		
		Matrix[] x = new Matrix[x0.length+1]; // the vertices of the simplex
		double[] fx = new double[x0.length+1]; // the values of f at the vertices
		for (int i = 0; i < x.length; i ++) {
			x[i] = new Matrix(new double[][] {x0.clone()}); // initialize the vertices as the guess
			if (i < x0.length)
				x[i].set(0, i, x[i].get(0, i)*1.9); // with multidimensional perturbations
			fx[i] = f.apply(x[i]); // and get the initial values
		}
		if (!Double.isFinite(fx[x0.length]))
			throw new IllegalArgumentException("Initial guess yielded bunk value");
		
		while (true) { // now for the iterative part
			int iWorst = NumericalMethods.argmax(fx);
			int iNext = NumericalMethods.argpenmax(fx);
			int iBest = NumericalMethods.argmin(fx);
			
			if (Math.abs((fx[iWorst] - fx[iBest])/fx[iBest]) < tol) // if these are all basically the same
				return x[iBest].values[0]; // terminate and take the best one we have
			
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
			Function<double[], Double> func, List<Function<double[], Double>> grad, double[] x0, double[] scale, double tol) {
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
			for (int i = 0; i < x0.length; i += n) { // iterate through the coordinates
				double[] y0 = Arrays.copyOfRange(x, i, i+n); // extract slices of the array
				final int I = i;
				double[] yOpt = minimizeNelderMead((y) -> { // then run Nelder-Mead on that slice
					double[] xp = x.clone();
					System.arraycopy(y, 0, xp, I, n);
					return func.apply(xp);
				}, y0, tol);
				System.arraycopy(yOpt, 0, x, i, n);
			} // after going through each coordinate
			double fxp = func.apply(x); // check how much change that whole cycle got us
			if (Math.abs((fxp - fx)/fx) < tol)
				return x;
			else
				fx = fxp;
		}
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
				H0 = 1;
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
			
			if (Math.abs((Ui - Uf)/Ui) < tol) { // STEP 4: stop condition
				return x.T().values[0]; // if the energy isn't really changing, then we're done
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
	public static class Matrix {
		
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
		 * Instantiate a matrix based on a 1D array given dimensions.
		 * @param n - The height
		 * @param m - The width
		 * @param values - The values from left to right then top to bottom
		 */
		public Matrix(int n, int m, double... values) {
			if (values.length != n*m)
				throw new IllegalArgumentException(values.length+" values do not fit in a "+n+"x"+m+" matrix.");
			this.values = new double[n][m];
			for (int i = 0; i < values.length; i ++)
				this.set(i/m, i%m, values[i]);
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
						prod.set(i,j,prod.get(i,j) + this.get(i, k)*that.get(k, j));
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
		 * Set a single scalar value.
		 * @param i
		 * @param j
		 * @param a
		 */
		public void set(int i, int j, double a) {
			this.values[i][j] = a;
		}
		
		public String toString() {
			String str = "Matrix("+this.getN()+", "+this.getM()+",\n  ";
			for (int i = 0; i < this.getN(); i ++) {
				for (int j = 0; j < this.getM(); j ++)
					str += String.format("%- 6.4f, ", this.get(i, j));
				str += "\n  ";
			}
			return str.substring(0, str.length()-5) + ")";
		}
	}
	
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		Function<double[], Double> simionescu = (v) -> {
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
		Function<double[], Double> himmelblau = (v) -> {
			double x = v[0], y = v[1];
			double z = Math.pow(x*x + y - 11, 2) + Math.pow(x + y*y - 7, 2) + 1;
			System.out.printf("[%.4f, %.4f, %.4f],\n", x, y, z);
			return z;
		};
		
//		System.out.println(Arrays.toString(minimizeNelderMead(
//				simionescu, new double[] {.5,.5}, 1e-8)));
		System.out.println(Arrays.toString(minimizeCoordinateDescent(
				simionescu, simionescuGradient, new double[] {.5, .5}, new double[] {1,1}, 1e-8)));
//		System.out.println(Arrays.toString(minimizeCoordinateDescent(
//				himmelblau, new double[] {0,0}, new double[] {1,1}, 1e-8)));
//		System.out.println(Arrays.toString(minimizeLBFGS(
//				himmelblau, new double[] {.5, .5}, 1e-8)));
	}


}
