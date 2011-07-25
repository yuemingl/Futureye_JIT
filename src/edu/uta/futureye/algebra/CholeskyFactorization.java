package edu.uta.futureye.algebra;

import edu.uta.futureye.algebra.intf.Matrix;

public class CholeskyFactorization {

	/**
	 * A=A' and x'Ax>0 (Positive definite)
	 * 
	 * @param A
	 * @param L
	 */
	public static void Cholesky(Matrix AA, Matrix L) {
		Matrix A = AA.copy();
		int N = A.getRowDim();
		
		for(int n=1; n<=N; n++) {
			double lnn = Math.sqrt(A.get(n, n));
			L.set(n, n, lnn);
			for(int i=n; i<=N; i++) {
				L.set(i, n, A.get(i, n)/lnn);
			}
			for(int i=n; i<=N; i++) {
				for(int j=n; j<=N; j++) {
					double lin = L.get(i, n);
					double ljn = L.get(j, n);
					A.add(i, j, -lin*ljn);
				}
			}
		}
	}
	
	public static void test1() {
		SparseMatrix A = new SparseMatrix(3,3);
		SparseMatrix L = new SparseMatrix(3,3);
		
		double[][] data = {{25,15,-5},{15,18,0},{-5,0,11}};
		for(int i=0;i<data.length;i++) {
			for(int j=0;j<data[i].length;j++)
				A.set(i+1, j+1, data[i][j]);
		}
		A.print();
		Cholesky(A,L);
		A.print();
		L.print();		
	}
	
	public static void test2() {
		SparseMatrix A = new SparseMatrix(6,6);
		SparseMatrix L = new SparseMatrix(6,6);
		
		double[][] data = 
	      {{1,     1,     1,     1,     1,     1},
	       {1,     2,     3,     4,     5,     6},
	       {1,     3,     6,    10,    15,    21},
	       {1,     4,    10,    20,    35,    56},
	       {1,     5,    15,    35,    70,   126},
	       {1,     6,    21,    56,   126,   252}};
	       
	       
		for(int i=0;i<data.length;i++) {
			for(int j=0;j<data[i].length;j++)
				A.set(i+1, j+1, data[i][j]);
		}
		A.print();
		Cholesky(A,L);
		A.print();
		L.print();		
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		//test1();
		test2();
		
	}

}
