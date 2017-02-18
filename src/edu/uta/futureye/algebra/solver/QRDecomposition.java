/**
 * Copyright (c) 2010, nkliuyueming@gmail.com. All rights reserved.
 * 
 * 
 */
package edu.uta.futureye.algebra.solver;

import edu.uta.futureye.algebra.SparseMatrixRowMajor;
import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.SparseMatrix;

public class QRDecomposition {
	
	private static double getColumnNorm(Matrix A, int col) {
		int N = A.getRowDim();
		double norm = 0.0;
		for(int n=1; n<=N; n++) {
			double an = A.get(n, col);
			norm += an*an;
		}
		return Math.sqrt(norm);
	}

	/**
	 * A=Q*R using the Gram-Schmidt process
	 * where Q'Q=I
	 *  
	 * @param A
	 * @param Q
	 * @param R
	 */
	public static void QR(Matrix A, Matrix Q, Matrix R) {
		Matrix AA = A.copy();
		int N = AA.getRowDim();
		for(int n=1; n<=N; n++) {
			double rnn = getColumnNorm(AA, n);
			R.set(n, n, rnn);
			for(int i=1;i<=N;i++) {
				Q.set(i, n, AA.get(i, n)/rnn);
			}
			for(int i=n+1;i<=N;i++) {
				double R12j = 0.0;
				for(int j=1;j<=N;j++) {
					R12j += Q.get(j, n)*AA.get(j, i);
				}
				R.set(n, i, R12j);
			}
			for(int i=n+1;i<=N;i++) {
				for(int j=1;j<=N;j++) {
					double Q2R22ij = Q.get(j, n)*R.get(n, i);
					AA.add(j, i, -Q2R22ij);
				}
			}	
		}
	}
	
	public static void test1() {
		SparseMatrix A = new SparseMatrixRowMajor(3,2);
		SparseMatrix Q = new SparseMatrixRowMajor(3,2);
		SparseMatrix R = new SparseMatrixRowMajor(2,2);
		
		double[][] data = {{3,-6},{4,-8},{0,1}};
		for(int i=0;i<data.length;i++) {
			for(int j=0;j<data[i].length;j++)
				A.set(i+1, j+1, data[i][j]);
		}
		A.print();
		QR(A,Q,R);
		A.print();
		Q.print();
		R.print();		
	}
	
	public static void test2() {
		SparseMatrix A = new SparseMatrixRowMajor(4,3);
		SparseMatrix Q = new SparseMatrixRowMajor(4,3);
		SparseMatrix R = new SparseMatrixRowMajor(3,3);
		
		double[][] data = {{9,0,26},{12,0,-7},{0,4,4},{0,-3,-3}};
		for(int i=0;i<data.length;i++) {
			for(int j=0;j<data[i].length;j++)
				A.set(i+1, j+1, data[i][j]);
		}
		A.print();
		QR(A,Q,R);
		A.print();
		Q.print();
		R.print();		
	}	
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		//test1();
		test2();
		
	}

}
