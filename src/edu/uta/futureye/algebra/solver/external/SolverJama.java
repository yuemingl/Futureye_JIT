/**
 * Copyright (c) 2010, nkliuyueming@gmail.com. All rights reserved.
 * 
 * 
 */
package edu.uta.futureye.algebra.solver.external;

import edu.uta.futureye.algebra.FullMatrix;
import edu.uta.futureye.algebra.SparseMatrixRowMajor;
import edu.uta.futureye.algebra.intf.SparseMatrix;
import edu.uta.futureye.algebra.solver.LUDecomposition;

/**
 * Jama Interface 
 * 
 * @author liuyueming
 *
 */
public class SolverJama {
	
	/**
	 * Solve <tt>A*X = B</tt> with LU Decomposition.
	 * @param A
	 * @param B
	 * @return X
	 */
	public FullMatrix solve(FullMatrix A,FullMatrix B) {
		
		Jama.Matrix JA = new Jama.Matrix(A.getData());
		Jama.Matrix JB = new Jama.Matrix(B.getData());
		
		long begin,end;
		begin = System.currentTimeMillis();
		Jama.LUDecomposition lu = JA.lu();
		end = System.currentTimeMillis();
		System.out.println("Jama LU="+(end-begin)+"ms");
		
		begin = System.currentTimeMillis();
		Jama.Matrix JX = lu.solve(JB);
		end = System.currentTimeMillis();
		System.out.println("Jama Solve="+(end-begin)+"ms");
		
		FullMatrix rlt = new FullMatrix(JX.getArray(),false);
		return rlt;
	}
	
	public static void test() {
		int N = 4024;
		int NE = 1;
		double[][] dataA = new double[N][N];
		double[][] dataB = new double[N][NE];
		double[] alpha = new double[N];
		for(int i=N-1;i>=0;i--) {
			alpha[i] = 0.99+(i*1.0)/N/1000.0;
		}
		for(int i=0;i<N;i++) {
			dataA[i][0] = 1.0;
		}
		for(int j=1;j<N;j++) {
			for(int i=0;i<N;i++) {
				dataA[i][j] = dataA[i][j-1]*alpha[i];
			}
		}
		for(int i=N-1;i>0;i--) {
			for(int j=NE-1;j>0;j--) {
				dataB[i][j] = i;
			}
		}
		FullMatrix A = new FullMatrix(dataA,true);
		FullMatrix B = new FullMatrix(dataB,true);
		SolverJama solver = new SolverJama();
		System.out.println("Begin Jama solver...");
		long begin = System.currentTimeMillis();
		FullMatrix X = solver.solve(A, B);
		long end = System.currentTimeMillis();
		System.out.println("Jama time used: "+(end-begin));

		FullMatrix fL = new FullMatrix(N,N);
		FullMatrix fU = new FullMatrix(N,N);
		FullMatrix fX = new FullMatrix(N,NE);
		FullMatrix fF = new FullMatrix(dataB,true);
		
		SparseMatrix P = new SparseMatrixRowMajor(N,N);
		System.out.println("Begin FuturEye solver...");
		begin = System.currentTimeMillis();
		//LUDecomposition.LU(A, fL, fU, P);
		LUDecomposition.solve(A, fL, fU, P, fX, fF);
		end = System.currentTimeMillis();
		System.out.println("FuturEye time used: "+(end-begin));
		
	}
	
	public static void main(String[] args) {
		test();
	}
}
