package edu.uta.futureye.algebra.solver.eigen;

import edu.uta.futureye.algebra.FullMatrix;
import edu.uta.futureye.algebra.FullVector;
import edu.uta.futureye.algebra.intf.AlgebraMatrix;
import edu.uta.futureye.algebra.intf.AlgebraVector;
import edu.uta.futureye.algebra.solver.Solver;
import edu.uta.futureye.util.Constant;
import edu.uta.futureye.util.FutureyeException;

public class PowerMethod {
	protected static int maxIter = 10000;
	protected static double eps = Constant.eps;
	public static boolean debug = true;
	
	/**
	 * 
	 * @param nmax
	 */
	public static void setMaximunIteration(int nmax) {
		maxIter = nmax;
	}
	
	/**
	 * 
	 * @param eps
	 */
	public static void setEpsilon(double eps) {
		PowerMethod.eps = eps;
	}	
	
	/**
	 * Compute largest eigenvalue in magnitude of A
	 * 
	 * @param A
	 * @param x ==null, just return eigenvalue
	 *        x !=null, x will contains the corresponding eigenvector
	 * @return
	 */
	public static double largestEigenvalue(AlgebraMatrix A, AlgebraVector x) {
		FullVector x0 = new FullVector(A.getColDim());
		FullVector x1 = new FullVector(A.getColDim());
		x0.setRandom(1.0, 0.0);
		double lmd0 = Double.MIN_VALUE;
		double lmd = 0;
		int i=maxIter;
		while( --i>0 ) {
			A.mult(x0, x1);
			/* *******************************************************
			 * Another way of finding eigenvalue is by Rayleigh quotient
			 * lmd = x0'*A*x0/x0'x0 = x0'*x1
			 **********************************************************/
			lmd = x1.get(1)/x0.get(1); //eigenvalue
			x1.scale(1.0/x1.normInf());
			FullVector tmp = x0;
			x0 = x1;
			x1 = tmp;
			if(Math.abs(lmd-lmd0) < eps) {
				break;
			}
			lmd0 = lmd;
		}
		if(x != null) x.set(x0);
		if(debug) {
			System.out.println("Number of iterations: "+(maxIter-i));
			//res = A*x0 - lmd*x0
			A.mult(x0, x1);
			FullVector res = x0.scale(lmd).axpy(-1.0, x1);
			System.out.println("Residual norm2(A*x-lmd*x) = "+res.norm2());
		}
		if(i == 0) {
			throw new FutureyeException("Maximun number of iterations reached: "+maxIter);
		}
		return lmd;
	}
	
	/**
	 * Compute smallest eigenvalue in magnitude of A
	 * 
	 * @param A
	 * @param x ==null, just return eigenvalue
	 *        x !=null, x will contains the corresponding eigenvector
	 * @return
	 */
	public static double smallestEigenvalue(AlgebraMatrix A, AlgebraVector x) {
		FullVector x0 = new FullVector(A.getColDim());
		FullVector x1 = new FullVector(A.getColDim());
		x0.setRandom(1.0, 0.0);
		double lmd0 = Double.MAX_VALUE;
		double lmd = 0;
		int i=maxIter;
		Solver sol = new Solver();
		sol.epsAbsIterMax = eps;
		sol.epsAbsIterMin = eps;
		while( --i>0 ) {
			sol.solveCGS(A, x0, x1);
			lmd = x1.get(1)/x0.get(1);
			x1.scale(1.0/x1.normInf());
			FullVector tmp = x0;
			x0 = x1;
			x1 = tmp;
			if(Math.abs(lmd-lmd0) < eps) {
				break;
			}
			lmd0 = lmd;
		}
		lmd = 1.0/lmd;
		if(x != null) x.set(x0);
		if(debug) {
			System.out.println("Number of iterations: "+(maxIter-i));
			//res = A*x0 - lmd*x0
			A.mult(x0, x1);
			FullVector res = x0.scale(lmd).axpy(-1.0, x1);
			System.out.println("Residual norm2(A*x-lmd*x) = "+res.norm2());
		}
		if(i == 0) {
			throw new FutureyeException("Maximun number of iterations reached: "+maxIter);
		}
		return lmd;
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		double[][] data1 = {{1,3,-7},{-3,4,1},{2,-5,3}};
		//   3.9893 + 5.5601i
		//   3.9893 - 5.5601i
		//   0.0214
		
		double[][] data2 = {{6,5,-5},{2,6,2},{2,5,-1}};
	    //  8.6235
	    //  4.0000
	    // -1.6235
		
		FullMatrix mat = new FullMatrix(data2,false);
		FullVector x = new FullVector(mat.getColDim());
		double lmd1 = largestEigenvalue(mat,x);
		System.out.println(lmd1);
		
		x.print();
		double lmd2 = smallestEigenvalue(mat, x);
		System.out.println(lmd2);
		x.print();

	}

}
