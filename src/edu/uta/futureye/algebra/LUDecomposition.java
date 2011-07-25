package edu.uta.futureye.algebra;

import edu.uta.futureye.algebra.intf.Matrix;

public class LUDecomposition {

	/*
	 * A=L*U
	 */
	public static void LU(Matrix A, Matrix L, Matrix U) {
		int N = A.getRowDim();
		
		for(int n=1; n<=N; n++) {
			for(int j=n; j<=N; j++) {
				double sum1 = 0.0;
				for(int k=1; k<=n-1; k++) {
					sum1 += L.get(n, k)*U.get(k, j);
				}
				double anj = A.get(n, j);
				U.set(n, j, anj-sum1);
			}
			for(int i=n+1; i<=N; i++) {
				double sum2 = 0.0;
				for(int k=1; k<=n-1; k++) {
					sum2 += L.get(i, k)*U.get(k, n);
				}
				double ain = A.get(i, n);
				L.set(i, n, (ain-sum2)/U.get(n, n));
			}
			L.set(n, n, 1.0);
		}
	}
	
	/*
	 * A=P*L*U
	 * 
	 */
	public static void LU(Matrix AA, Matrix L, Matrix U, Matrix P) {
		Matrix A = AA.copy();
		int N = A.getRowDim();
		int[] VP = new int[N+1];
		for(int n=1; n<=N; n++) {
			VP[n] = n;
		}
		for(int n=1; n<=N; n++) {
			//判断Unn的值是否为0，为0的话进行行对换，直到找到Unn不为0的行
			for(int nn=n; nn<=N; nn++) {
				double sum1 = 0.0;
				for(int k=1; k<=n-1; k++) {
					sum1 += L.get(nn, k)*U.get(k, nn);
				}
				double ann = A.get(nn, n);
				if(Math.abs(ann-sum1) > 1e-12) {
					if(nn!=n) {
						int tmp = VP[nn];
						VP[nn] = VP[n];
						VP[n] = tmp;
						double vtmp;
						for(int c=1;c<=N;c++) {
							vtmp = A.get(n, c);
							A.set(n, c, A.get(nn, c));
							A.set(nn, c, vtmp);
						}
						for(int c=1;c<n;c++) {
							vtmp = L.get(n, c);
							L.set(n, c, L.get(nn, c));
							L.set(nn, c, vtmp);
						}
					}
					break;
				}
			}
			for(int j=n; j<=N; j++) {
				double sum1 = 0.0;
				for(int k=1; k<=n-1; k++) {
					sum1 += L.get(n, k)*U.get(k, j);
				}
				double anj = A.get(n, j);
				U.set(n, j, anj-sum1);
			}
			for(int i=n+1; i<=N; i++) {
				double sum2 = 0.0;
				for(int k=1; k<=n-1; k++) {
					sum2 += L.get(i, k)*U.get(k, n);
				}
				double ain = A.get(i, n);
				L.set(i, n, (ain-sum2)/U.get(n, n));
			}
			L.set(n, n, 1.0);
		}
		
		//P
		for(int i=1;i<=N;i++) {
			P.set(i, VP[i], 1.0);
		}
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		SparseMatrix A = new SparseMatrix(3,3);
		SparseMatrix L = new SparseMatrix(3,3);
		SparseMatrix U = new SparseMatrix(3,3);
		SparseMatrix P = new SparseMatrix(3,3);
		
		//double[][] data = {{8,2,9},{4,9,4},{6,7,9}};
		double[][] data = {{1,0,0},{0,0,2},{0,1,-1}};
		for(int i=0;i<data.length;i++) {
			for(int j=0;j<data[i].length;j++)
				A.set(i+1, j+1, data[i][j]);
		}
		A.print();
		LU(A,L,U,P);
		A.print();
		L.print();
		U.print();
		P.print();


	}

}
