/**
 * Copyright (c) 2010, nkliuyueming@gmail.com. All rights reserved.
 * 
 * 
 */
package edu.uta.futureye.algebra.solver;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import edu.uta.futureye.algebra.FullMatrix;
import edu.uta.futureye.algebra.FullVector;
import edu.uta.futureye.algebra.SparseMatrixColMajor;
import edu.uta.futureye.algebra.SparseMatrixRowMajor;
import edu.uta.futureye.algebra.SparseVectorHashMap;
import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.SparseMatrix;
import edu.uta.futureye.algebra.intf.SparseVector;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.util.FutureyeException;

public class LUDecomposition {
	public static double eps = 1E-50;

	/**
	 * A=L*U  without permutation of L (this function is just for test only)
	 * 
	 * @param A (Input)
	 * @param L (Output) Lower triangular matrix with ones on its diagonal
	 * @param U (Output) Upper triangular matrix
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
	
	/**
	 * LU Decomposition without pivot
	 * 
	 * @param A
	 * @param L
	 * @param U
	 */
	public static void LU(FullMatrix A, FullMatrix L, FullMatrix U) {
		int N = A.getRowDim();
		double[][] dA = A.getData();
		double[][] dL = L.getData();
		double[][] dU = U.getData();
		
		for(int n=0; n<=N-1; n++) {
			for(int j=n; j<N; j++) {
				double sum1 = 0.0;
				for(int k=0; k<=n-1; k++) {
					sum1 += dL[n][k]*dU[k][j];
				}
				double anj = dA[n][j];
				dU[n][j] = anj-sum1;
			}
			for(int i=n+1; i<N; i++) {
				double sum2 = 0.0;
				for(int k=0; k<=n-1; k++) {
					sum2 += dL[i][k]*dU[k][n];
				}
				double ain = dA[i][n];
				dL[i][n] = (ain-sum2)/dU[n][n];
			}
			dL[n][n] = 1.0;
		}
	}
	
	/**
	 * A=P*L*U
	 * 
	 * @param A
	 * @param L (Output) Lower triangular matrix with ones on its diagonal
	 * @param U (Output) Upper triangular matrix
	 * @param P (Output) Permutation matrix
	 */
	public static void LU(Matrix A, Matrix L, Matrix U, Matrix P) {
		Matrix AA = A.copy();
		int N = AA.getRowDim();
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
				double ann = AA.get(nn, n);
				if(Math.abs(ann-sum1) > 1e-12) {
					if(nn!=n) {
						int tmp = VP[nn];
						VP[nn] = VP[n];
						VP[n] = tmp;
						double vtmp;
						for(int c=1;c<=N;c++) {
							vtmp = AA.get(n, c);
							AA.set(n, c, AA.get(nn, c));
							AA.set(nn, c, vtmp);
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
				double anj = AA.get(n, j);
				U.set(n, j, anj-sum1);
			}
			for(int i=n+1; i<=N; i++) {
				double sum2 = 0.0;
				for(int k=1; k<=n-1; k++) {
					sum2 += L.get(i, k)*U.get(k, n);
				}
				double ain = AA.get(i, n);
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
	 * A=P*L*U
	 * Sparse LU decomposition (dev version)
	 * 
	 * @param A (Input) Sparse coefficient matrix
	 * @param L (Output) Lower triangular matrix with ones on its diagonal
	 * @param U (Output) Upper triangular matrix
	 * @param P (Output) Permutation matrix
	 */
	public static void LU2(SparseMatrixRowMajor A, SparseMatrixRowMajor L, SparseMatrix U, SparseMatrix P) {
		SparseMatrixRowMajor AA = A.copy();
		SparseMatrixColMajor UU = new SparseMatrixColMajor(U.getRowDim(),U.getColDim());
		
		int N = AA.getRowDim();
		int[] VP = new int[N+1];
		//for(int n=1; n<=N; n++) {
		for(int n=N+1; --n>=1;) {
			VP[n] = n;
		}
		
		long time1=0,time2=0,time3=0;
		long begin,end;
		int swap = 0;
		double[] cache = new double[N];
		for(int i=N; --i>=0;) cache[i] = 0.0;
		
		for(int n=1; n<=N; n++) {
			//判断Unn的值是否为0，为0的话进行行对换，直到找到Unn不为0的行
			begin = System.currentTimeMillis();
			for(int nn=n; nn<=N; nn++) {
				double sum1 = L.mult(UU, nn, nn);
				double ann = AA.get(nn, n);
				if(Math.abs(ann-sum1) > 1e-50) {
					if(nn!=n) {
						int tmp = VP[nn];
						VP[nn] = VP[n];
						VP[n] = tmp;
						AA.swapRow(n, nn);
						L.swapRow(n, nn);
						swap++;
					}
					break;
				}
			}
			end = System.currentTimeMillis();
			time1 += end-begin;
			
			begin = System.currentTimeMillis();
			if(n>1) {
				Map<Integer, Double> LnLast = L.getAll().get(n-1);
				if(LnLast != null) {
				for(Entry<Integer,Double> e : LnLast.entrySet())
					cache[e.getKey()-1] = 0.0;
				}
			}
			Map<Integer, Double> Ln = L.getAll().get(n);
			if(Ln != null) {
			for(Entry<Integer,Double> e : Ln.entrySet())
				cache[e.getKey()-1] = e.getValue();
			}
			
			for(int j=n; j<=N; j++) {
				//double sum1 = L.mult(UU, n, j);
				double sum1 = UU.multColumn(cache, j);
				double anj = AA.get(n, j);
				UU.set(n, j, anj-sum1);
			}
			end = System.currentTimeMillis();
			time2 += end-begin;
			
			begin = System.currentTimeMillis();
			double Unn = UU.get(n, n);
			for(int i=n+1; i<=N; i++) {
//				if(L.getAll().get(i)!=null && UU.getAll().get(n)!=null) {
//					System.out.println("L="+L.getAll().get(i).size()+" U="+UU.getAll().get(n).size());
//				}
				double sum2 = L.mult(UU, i, n);
				double ain = AA.get(i, n);
				L.set(i, n, (ain-sum2)/Unn);
			}
			end = System.currentTimeMillis();
			time3 += end-begin;
			
			L.set(n, n, 1.0);
		}
		
		System.out.println("Sparse LU time1="+time1 + " swap="+swap);
		System.out.println("Sparse LU time2="+time2);
		System.out.println("Sparse LU time3="+time3);
		
		//P
		for(int i=1;i<=N;i++) {
			P.set(i, VP[i], 1.0);
		}
		
		//U
		Map<Integer, Map<Integer, Double>> dataU = UU.getAll();
		for(Entry<Integer, Map<Integer, Double>> e1 : dataU.entrySet()) {
			Map<Integer, Double> col = e1.getValue();
			int nCol = e1.getKey();
			for(Entry<Integer, Double> e2 : col.entrySet()) {
				int nRow = e2.getKey();
				U.set(nRow, nCol, e2.getValue());
			}
		}
		
	}
	
	/**
	 * A=P*L*U
	 * Sparse LU decomposition
	 * 
	 * @param A (Input) Sparse coefficient matrix
	 * @param L (Output) Lower triangular matrix with ones on its diagonal
	 * @param U (Output) Upper triangular matrix
	 * @param P (Output) Permutation matrix
	 */
	public static void LU(SparseMatrixRowMajor A, SparseMatrixRowMajor L, SparseMatrix U, SparseMatrix P) {
		SparseMatrixRowMajor AA = A.copy();
		SparseMatrixColMajor UU = new SparseMatrixColMajor(U.getRowDim(),U.getColDim());
		
		int N = AA.getRowDim();
		int[] VP = new int[N+1];
		//for(int n=1; n<=N; n++) {
		for(int n=N+1; --n>=1;) {
			VP[n] = n;
		}
		
		long time1=0,time2=0,time3=0;
		long begin,end;
		int swap = 0;
		double[] cache = new double[N];
		for(int i=N; --i>=0;) cache[i] = 0.0;
		List<Integer> nonzeroIndex = new ArrayList<Integer>();
		
		for(int n=1; n<=N; n++) {
			//判断Unn的值是否为0，为0的话进行行对换，直到找到Unn不为0的行
			begin = System.currentTimeMillis();
			for(int nn=n; nn<=N; nn++) {
				double sum1 = L.mult(UU, nn, nn);
				double ann = AA.get(nn, n);
				if(Math.abs(ann-sum1) > 1e-50) {
					if(nn!=n) {
						int tmp = VP[nn];
						VP[nn] = VP[n];
						VP[n] = tmp;
						AA.swapRow(n, nn);
						L.swapRow(n, nn);
						swap++;
					}
					break;
				}
			}
			end = System.currentTimeMillis();
			time1 += end-begin;
			
			begin = System.currentTimeMillis();
			nonzeroIndex.clear();
			for(int j=1; j<n; j++) {
				double sum1 = UU.multColumn(cache,nonzeroIndex, j);
				double Lnj = (AA.get(n, j)-sum1)/UU.get(j, j);
				if(Math.abs(Lnj) > Matrix.zeroEps) {
					nonzeroIndex.add(j);
					cache[j-1] = Lnj;
					L.set(n, j, Lnj);
				}
			}
			end = System.currentTimeMillis();
			time3 += end-begin;
			
			begin = System.currentTimeMillis();
			for(int j=n; j<=N; j++) {
				double sum1 = UU.multColumn(cache,nonzeroIndex, j);
				double anj = AA.get(n, j);
				UU.set(n, j, anj-sum1);
			}
			end = System.currentTimeMillis();
			time2 += end-begin;
			
			L.set(n, n, 1.0);
		}
		
		System.out.println("FuturEye Sparse LU time1/3="+time1 + " swap="+swap);
		System.out.println("FuturEye Sparse LU time2/3="+time2);
		System.out.println("FuturEye Sparse LU time3/3="+time3);
		
		//P
		for(int i=1;i<=N;i++) {
			P.set(i, VP[i], 1.0);
		}
		
		//U
		Map<Integer, Map<Integer, Double>> dataU = UU.getAll();
		for(Entry<Integer, Map<Integer, Double>> e1 : dataU.entrySet()) {
			Map<Integer, Double> col = e1.getValue();
			int nCol = e1.getKey();
			for(Entry<Integer, Double> e2 : col.entrySet()) {
				int nRow = e2.getKey();
				U.set(nRow, nCol, e2.getValue());
			}
		}
		
	}
	
	/**
	 * A=P*L*U
	 * Full matrix LU decomposition with pivot
	 * 
	 * @param A (Input) Full coefficient matrix
	 * @param L (Output) Lower triangular matrix with ones on its diagonal
	 * @param U (Output) Upper triangular matrix
	 * @param P (Output) Permutation matrix
	 */
	public static void LU(FullMatrix A, FullMatrix L, FullMatrix U, SparseMatrix P) {
		FullMatrix AA = A.copy();
		int N = AA.getRowDim();
		int[] VP = new int[N];
		for(int n=0; n<N; n++) {
			VP[n] = n;
		}
		double[][] dA = AA.getData();
		double[][] dL = L.getData();
		double[][] dU = U.getData();
		double[][] dUT = new double[N][N];
		for(int i=0;i<N;i++) {
			for(int j=0;j<N;j++) {
				dUT[j][i] = dU[i][j];
			}
		}
		
		//long begin=0,end=0;
		//long time1=0,time2=0,time3=0;
		for(int n=0; n<N; n++) {
			//判断Unn的值是否为0，为0的话进行“行对换”，直到找到Unn不为0的行
			boolean bFind = false;
			for(int nn=n; nn<N; nn++) {
				double sum1 = 0.0;
				double[] _dUTnn = dUT[nn];
				//begin = System.currentTimeMillis();
				double[] _dLnn = dL[nn];
				for(int k=0; k<=n-1; k++) {
					sum1 += _dLnn[k]*_dUTnn[k];
				}
				//end = System.currentTimeMillis();
				//time1 += end-begin;
				
				//begin = System.currentTimeMillis();
				double ann = dA[nn][n];
				if(Math.abs(ann-sum1) > eps) {
					if(nn!=n) {
						int tmp = VP[nn];
						VP[nn] = VP[n];
						VP[n] = tmp;
						double vtmp;
						for(int c=0;c<N;c++) {
							vtmp = dA[n][c];
							dA[n][c] = dA[nn][c];
							dA[nn][c] = vtmp;
						}
						for(int c=0;c<n;c++) {
							vtmp = dL[n][c];
							dL[n][c] = dL[nn][c];
							dL[nn][c] = vtmp;
						}
					}
					//System.out.println(nn);
					bFind = true;
					break;
				}
				
				//end = System.currentTimeMillis();
				//time3 += end-begin;
			}
			if(!bFind)
				throw new FutureyeException("Matrix is singular!");
			
			double[] _dAn = dA[n];
			double[] _dLn = dL[n];
			//begin = System.currentTimeMillis();
			for(int j=n; j<N; j++) {
				double sum1 = 0.0;
				double[] _dUTj = dUT[j];
				for(int k=0; k<=n-1; k++) {
					sum1 += _dLn[k]*_dUTj[k];
				}
				double anj = _dAn[j];
				dUT[j][n] = anj-sum1;
			}
			for(int i=n+1; i<N; i++) {
				double sum2 = 0.0;
				double[] _dLi = dL[i];
				double[] _dUTn = dUT[n];
				for(int k=0; k<=n-1; k++) {
					sum2 += _dLi[k]*_dUTn[k];
				}
				double ain = dA[i][n];
				dL[i][n] = (ain-sum2)/dUT[n][n];
			}
			//end = System.currentTimeMillis();
			//time2 += end-begin;
			dL[n][n] = 1.0;
		}
		//System.out.println("Time1="+time1);
		//System.out.println("Time2="+time2);
		//System.out.println("Time3="+time3);
		
		for(int i=0;i<N;i++) {
			for(int j=0;j<N;j++) {
				dU[i][j] = dUT[j][i];
			}
		}
		//P
		for(int i=1;i<=N;i++) {
			P.set(i, VP[i-1]+1, 1.0);
		}
	}
	
	/**
	 * Solve U*x=f
	 * 
	 * @param U (Input) Upper triangular matrix
	 * @param x (Input) Unknown vector 
	 * @param f (Input) Right hand side
	 */
	public static Vector solveUx(Matrix U, Vector x, Vector f) {
		int n = U.getRowDim();
		for(int i=n;i>0;i--) {
			double sum = 0.0;
			for(int j=n;j>i;j--) {
				sum += U.get(i, j)*x.get(j);
			}
			double xi = (f.get(i) - sum)/U.get(i, i);
			x.set(i, xi);
		}
		return x;
	}
	
	/**
	 * Solve U*x=f
	 * 
	 * @param U (Input) Upper triangular matrix
	 * @param x (Input) Unknown vector 
	 * @param f (Input) Right hand side
	 */
	public static FullVector solveUx(FullMatrix U, FullVector x, FullVector f) {
		int N = U.getRowDim();
		double[][] dU = U.getData();
		double[] dx = x.getData();
		double[] df = f.getData();
		for(int i=N-1;i>=0;i--) {
			double sum = 0.0;
			double[] _dUi = dU[i];
			for(int j=N-1;j>i;j--) {
				sum += _dUi[j]*dx[j];
			}
			double xi = (df[i] - sum)/dU[i][i];
			dx[i] = xi;
		}
		return x;
	}
	
	/**
	 * Solve U*x=f
	 * 
	 * @param U (Input) Upper triangular matrix
	 * @param x (Input) Unknown vector 
	 * @param f (Input) Right hand side
	 */
	public static FullMatrix solveUx(FullMatrix U, FullMatrix X, FullMatrix F) {
		int N = U.getRowDim();
		double[][] dU = U.getData();
		double[][] dX = X.getData();
		double[][] dF = F.getData();
		int nFCol = F.getColDim();
		double[] sum = new double[nFCol];
		for(int i=N-1;i>=0;i--) {
			for(int s=nFCol-1;s>=0;s--) sum[s] = 0.0;
			double[] _dUi = dU[i];
			for(int j=N-1;j>i;j--) {
				double _dUij = _dUi[j];
				double[] _dXj = dX[j];
				for(int k=0;k<nFCol;k++) {//循环--右端项列数
					sum[k] += _dUij*_dXj[k];
				}
			}
			double[] _dFi = dF[i];
			double[] _dXi = dX[i];
			for(int k=0;k<nFCol;k++) {
				_dXi[k] = (_dFi[k] - sum[k])/dU[i][i];
			}
		}
		return X;
	}	
	
	/**
	 * Solve U*x=f
	 * 
	 * @param U (Input) Upper triangular matrix
	 * @param x (Input) Unknown vector 
	 * @param f (Input) Right hand side
	 */
	public static SparseMatrix solveUx(SparseMatrixRowMajor U, SparseMatrix X, SparseMatrix F) {
		int N = U.getRowDim();
		int nFCol = F.getColDim();
		SparseMatrixColMajor XX = new SparseMatrixColMajor(X.getRowDim(),X.getColDim());
		for(int i=N+1;--i>0;) {
			for(int k=1;k<=nFCol;k++) {//循环--右端项列数
				double sum = U.mult(XX, i,k);
				XX.set(i, k, (F.get(i, k)-sum)/U.get(i, i));
			}
		}
		//X
		Map<Integer, Map<Integer, Double>> dataU = XX.getAll();
		for(Entry<Integer, Map<Integer, Double>> e1 : dataU.entrySet()) {
			Map<Integer, Double> col = e1.getValue();
			int nCol = e1.getKey();
			for(Entry<Integer, Double> e2 : col.entrySet()) {
				int nRow = e2.getKey();
				X.set(nRow, nCol, e2.getValue());
			}
		}
		return X;
	}	

	/**
	 * Solve U*x=f
	 * 
	 * @param U (Input) Upper triangular matrix
	 * @param x (Input) Unknown vector 
	 * @param f (Input) Right hand side
	 */
	public static FullMatrix solveUx(SparseMatrixRowMajor U, FullMatrix X, FullMatrix F) {
		int N = U.getRowDim();
		int nFCol = F.getColDim();
		double[][] dX = X.getData();
		double[][] dF = F.getData();
		double[] row = new double[N];
		int[] index = new int[N];
		double[] sums = new double[nFCol];
		int size;
		for(int i=N;--i>=0;) {
			size = 0;
			Map<Integer,Double> dMap = U.viewRow(i+1).getAll();
			if(dMap != null) {
				for(Entry<Integer,Double> e : dMap.entrySet()) {
					int idx = e.getKey()-1;
					if(idx!=i) {
						index[size] = idx;
						row[size] = e.getValue();
						size++;
					}
				}
			}
			double Uii = U.get(i+1, i+1);
			
			if(size>0) {
				for(int k=0;k<nFCol;k++) //循环--右端项列数
					sums[k] = row[0]*dX[index[0]][k];
				for(int s=1;s<size;s++)
					for(int k=0;k<nFCol;k++) //循环--右端项列数
						sums[k] += row[s]*dX[index[s]][k];
				for(int k=0;k<nFCol;k++) {//循环--右端项列数
					dX[i][k] = (dF[i][k]-sums[k])/Uii;
				}
			} else {
				for(int k=0;k<nFCol;k++)
					dX[i][k] = dF[i][k]/Uii;
			}
			
			
//			for(int k=0;k<nFCol;k++) {//循环--右端项列数
//				double sum = 0.0;
//				for(int s=0;s<size;s++)
//					sum += row[s]*dX[index[s]][k];
//				dX[i][k] = (dF[i][k]-sum)/Uii;
//			}
		}
		return X;
	}	
	
	/**
	 * Solve L*x = f
	 * 
	 * @param L (Input) Lower triangular matrix with ones on its diagonal
	 * @param x (Input) Unknown vector 
	 * @param f (Input) Right hand side
	 * @return
	 */
	public static Vector solveLx(Matrix L, Vector x, Vector f) {
		int N = L.getRowDim();
		for(int i=1;i<=N;i++) {
			double sum = 0.0;
			for(int j=1;j<i;j++) {
				sum += L.get(i, j)*x.get(j);
			}
			double xi = f.get(i) - sum;
			x.set(i, xi);
		}
		return x;		
	}

	/**
	 * Solve L*x = f
	 * 
	 * @param L (Input) Lower triangular matrix with ones on its diagonal
	 * @param x (Input) Unknown vector 
	 * @param f (Input) Right hand side
	 * @return
	 */
	public static FullVector solveLx(FullMatrix L, FullVector x, FullVector f) {
		int N = L.getRowDim();
		double[][] dL = L.getData();
		double[] dx = x.getData();
		double[] df = f.getData();
		for(int i=0;i<N;i++) {
			double sum = 0.0;
			double[] _dLi = dL[i];
			for(int j=0;j<i;j++) {
				sum += _dLi[j]*dx[j];
			}
			double xi = df[i] - sum;
			dx[i] = xi;
		}
		return x;
	}
	
	/**
	 * Solve L*x = f
	 * 
	 * @param L (Input) Lower triangular matrix with ones on its diagonal
	 * @param x (Input) Unknown vector 
	 * @param f (Input) Right hand side
	 * @return
	 */
	public static SparseMatrix solveLx(SparseMatrixRowMajor L, SparseMatrix X, SparseMatrix F) {
		int N = L.getRowDim();
		int nFCol = F.getColDim();
		SparseMatrixColMajor XX = new SparseMatrixColMajor(X.getRowDim(),X.getColDim());
//		long time1=0,time2=0;
//		long begin,end;
		
//		Map<Integer, Map<Integer, Double>> dL= L.getAll();
//		for(Entry<Integer, Map<Integer, Double>> e:dL.entrySet())
//			System.out.println(e.getValue().size());
		
		for(int i=1;i<=N;i++) {
			for(int k=1;k<=nFCol;k++) {//循环--右端项列数
				//begin = System.currentTimeMillis();
				double sum = L.mult(XX, i, k);
				//end = System.currentTimeMillis();
				//time1+=end-begin;
				
				//begin = System.currentTimeMillis();
				double v = F.get(i, k)-sum;
				XX.set(i, k, v);
				//end = System.currentTimeMillis();
				//time2+=end-begin;
			}
		}
		//System.out.println("solveLx s1="+(end-begin)+"ms");
		//System.out.println("solveLx s11="+time1+"ms");
		//System.out.println("solveLx s12="+time2+"ms");
		
		//X
		Map<Integer, Map<Integer, Double>> dataU = XX.getAll();
		for(Entry<Integer, Map<Integer, Double>> e1 : dataU.entrySet()) {
			Map<Integer, Double> col = e1.getValue();
			int nCol = e1.getKey();
			for(Entry<Integer, Double> e2 : col.entrySet()) {
				int nRow = e2.getKey();
				X.set(nRow, nCol, e2.getValue());
			}
		}
		return X;	
	}
	
	/**
	 * Solve L*x = f
	 * 
	 * @param L (Input) Lower triangular matrix with ones on its diagonal
	 * @param x (Input) Unknown vector 
	 * @param f (Input) Right hand side
	 * @return
	 */
	public static FullMatrix solveLx(FullMatrix L, FullMatrix X, FullMatrix F) {
		int N = L.getRowDim();
		double[][] dL = L.getData();
		double[][] dX = X.getData();
		double[][] dF = F.getData();
		int nFCol = F.getColDim();
		double[] sum = new double[nFCol];
		for(int i=0;i<N;i++) {
			for(int s=nFCol-1;s>=0;s--) sum[s] = 0.0;
			double[] _dLi = dL[i];
			for(int j=0;j<i;j++) {
				double[] _dXj = dX[j];
				double _dLij = _dLi[j];
				for(int k=0;k<nFCol;k++) {//循环--右端项列数
					sum[k] += _dLij*_dXj[k];
				}
			}
			double[] _dFi = dF[i];
			double[] _dXi = dX[i];
			for(int k=0;k<nFCol;k++) {
				_dXi[k] = _dFi[k] - sum[k];
			}
		}
		return X;	
	}
	
	/**
	 * Solve L*x = f
	 * 最快的
	 * 
	 * @param L (Input) Lower triangular matrix with ones on its diagonal
	 * @param x (Input) Unknown vector 
	 * @param f (Input) Right hand side
	 * @return
	 */
	public static FullMatrix solveLx(SparseMatrixRowMajor L, FullMatrix X, FullMatrix F) {
		int N = L.getRowDim();
		double[][] dX = X.getData();
		double[][] dF = F.getData();
		int nFCol = F.getColDim();
		double[] row = new double[N];
		int[] index = new int[N];
		double[] sums = new double[nFCol];
		int size;
		for(int i=0;i<N;i++) {
			Map<Integer,Double> dMap = L.viewRow(i+1).getAll();
			size = 0;
			for(Entry<Integer,Double> e : dMap.entrySet()) {
				double idx = e.getKey()-1;
				if(idx!=i) {
					index[size] = e.getKey()-1;
					row[size] = e.getValue();
					size++;
				}
			}
			if(size>0) {
				for(int k=0;k<nFCol;k++) //循环--右端项列数
					sums[k] = row[0]*dX[index[0]][k];
				for(int s=1;s<size;s++)
					for(int k=0;k<nFCol;k++) //循环--右端项列数
						sums[k] += row[s]*dX[index[s]][k];
				for(int k=0;k<nFCol;k++)
					dX[i][k] = dF[i][k]-sums[k];
			} else {
				for(int k=0;k<nFCol;k++) {//循环--右端项列数
					dX[i][k] = dF[i][k];
				}				
			}
			
//			for(int k=0;k<nFCol;k++) {//循环--右端项列数
//				double sum = 0.0;
//				for(int s=0;s<size;s++)
//					sum += row[s]*dX[index[s]][k];
//				dX[i][k] = dF[i][k]-sum;
//			}
		}
		return X;	
	}
	
	/**
	 * Solve P*x = f
	 * 
	 * @param P (Output) Permutation matrix
	 * @param x (Input) Unknown vector 
	 * @param f (Input) Right hand side
	 * @return
	 */
	public static Vector solvePx(SparseMatrix P, Vector x, Vector f) {
		Map<Integer, Map<Integer, Double>> p = P.getAll();
		int n = P.getRowDim();
		for(int i=1;i<=n;i++) {
			Map<Integer, Double> row = p.get(i);
			for(Entry<Integer, Double> c : row.entrySet()) {
				x.set(c.getKey(),f.get(i));
			}
		}
		return x;
	}
	
	/**
	 * Solve P*x = f
	 * 
	 * @param P (Output) Permutation matrix
	 * @param x (Input) Unknown vector 
	 * @param f (Input) Right hand side
	 * @return
	 */
	public static FullVector solvePx(SparseMatrix P, FullVector x, FullVector f) {
		Map<Integer, Map<Integer, Double>> p = P.getAll();
		int n = P.getRowDim();
		double[] dx = x.getData();
		double[] df = f.getData();
		for(int i=1;i<=n;i++) {
			Map<Integer, Double> row = p.get(i);
			for(Entry<Integer, Double> c : row.entrySet()) {
				dx[c.getKey()-1] = df[i-1];
			}
		}
		return x;
	}
	
	/**
	 * Solve P*x = f
	 * 
	 * @param P (Output) Permutation matrix
	 * @param x (Input) Unknown vector 
	 * @param f (Input) Right hand side
	 * @return
	 */
	public static FullMatrix solvePx(SparseMatrix P, FullMatrix X, FullMatrix F) {
		Map<Integer, Map<Integer, Double>> p = P.getAll();
		int n = P.getRowDim();
		int nFCol = F.getColDim();
		double[][] dX = X.getData();
		double[][] dF = F.getData();
		for(int i=1;i<=n;i++) {
			Map<Integer, Double> row = p.get(i);
			for(Entry<Integer, Double> c : row.entrySet()) {
				for(int k=0;k<nFCol;k++) {
					dX[c.getKey()-1][k] = dF[i-1][k];
				}
			}
		}
		return X;
	}
	
	/**
	 * Solve P*x = f
	 * 
	 * @param P (Output) Permutation matrix
	 * @param x (Input) Unknown vector 
	 * @param f (Input) Right hand side
	 * @return
	 */
	public static SparseMatrix solvePx(SparseMatrix P, SparseMatrix X, SparseMatrix F) {
		Map<Integer, Map<Integer, Double>> p = P.getAll();
		int N = P.getRowDim();
		int nFCol = F.getColDim();
		for(int i=1;i<=N;i++) {
			Map<Integer, Double> row = p.get(i);
			for(Entry<Integer, Double> c : row.entrySet()) {
				//for(int k=1;k<=nFCol;k++) {
				for(int k=nFCol+1; --k>=1;) {
					X.set(c.getKey(), k, F.get(i, k));
				}
			}
		}
		return X;
	}
	
	/**
	 * Solve A*X=F with LU decomposition, 
	 * where A can be any class that implements interface <tt>Matrix</tt>
	 * 
	 * @param A (Input) Coefficient matrix
	 * @param L (Output) Lower triangular matrix with ones on its diagonal
	 * @param U (Output) Upper triangular matrix
	 * @param P (Output) Permutation matrix
	 * @param x (Input) Unknown vector 
	 * @param f (Input) Right hand side
	 * @return
	 */
	public static Vector solve(Matrix A, Matrix L, Matrix U, SparseMatrix P,
			Vector x, Vector f) {
		LU(A, L, U, P);
		solvePx(P,x,f);
		Vector x2 = x.copy();
		solveLx(L,x2,x);
		solveUx(U,x,x2);
		return x;
	}
	
	/**
	 * Vector RHS
	 * Solve A*x=f with LU decomposition, where A is full matrix
	 * 
	 * @param A (Input) Coefficient matrix
	 * @param L (Output) Lower triangular matrix with ones on its diagonal
	 * @param U (Output) Upper triangular matrix
	 * @param P (Output) Permutation matrix
	 * @param x (Input) Unknown vector 
	 * @param f (Input) Right hand side
	 * @return
	 */
	public static FullVector solve(FullMatrix A, FullMatrix L, FullMatrix U, SparseMatrix P,
			FullVector x, FullVector f) {
		long begin,end;
		begin = System.currentTimeMillis();
		LU(A, L, U, P);
		end = System.currentTimeMillis();
		System.out.println("FuturEye Full LU="+(end-begin)+"ms");
		
		begin = System.currentTimeMillis();
		solvePx(P,x,f);
		end = System.currentTimeMillis();
		System.out.println("FuturEye Full SolvePx="+(end-begin));
		FullVector x2 = x.copy();
		
		begin = System.currentTimeMillis();
		solveLx(L,x2,x);
		end = System.currentTimeMillis();
		System.out.println("FuturEye Full SolveLx="+(end-begin));
		
		begin = System.currentTimeMillis();
		solveUx(U,x,x2);
		end = System.currentTimeMillis();
		System.out.println("FuturEye Full SolveUx="+(end-begin));
		return x;
	}
	
	/**
	 * Matrix RHS
	 * Solve A*X=F with LU decomposition, where A is full matrix
	 * 
	 * @param A (Input) Coefficient matrix
	 * @param L (Output) Lower triangular matrix with ones on its diagonal
	 * @param U (Output) Upper triangular matrix
	 * @param P (Output) Permutation matrix
	 * @param X (Input) Unknown vectors in matrix form
	 * @param F (Input) Right hand sides in matrix form
	 * @return
	 */
	public static FullMatrix solve(FullMatrix A, 
			FullMatrix L, FullMatrix U, SparseMatrix P,
			FullMatrix X, FullMatrix F) {
		long begin,end;
		begin = System.currentTimeMillis();
		LU(A, L, U, P);
		end = System.currentTimeMillis();
		System.out.println("FuturEye Full LU="+(end-begin)+"ms");
		
		begin = System.currentTimeMillis();
		solvePx(P,X,F);
		end = System.currentTimeMillis();
		System.out.println("FuturEye Full SolvePx="+(end-begin)+"ms");
		
		FullMatrix X2 = X.copy();
		begin = System.currentTimeMillis();
		solveLx(L,X2,X);
		end = System.currentTimeMillis();
		System.out.println("FuturEye Full SolveLx="+(end-begin)+"ms");
		
		begin = System.currentTimeMillis();
		solveUx(U,X,X2);
		end = System.currentTimeMillis();
		System.out.println("FuturEye Full SolveUx="+(end-begin)+"ms");
		
		return X;
	}
	
	/**
	 * Solve A*X=F with LU decomposition, where A is sparse matrix
	 * 
	 * @param A (Input) Coefficient matrix
	 * @param L (Output) Lower triangular matrix with ones on its diagonal
	 * @param U (Output) Upper triangular matrix
	 * @param P (Output) Permutation matrix
	 * @param X (Input) Unknown vectors in matrix form (column based)
	 * @param F (Input) Right hand sides in matrix form (column based)
	 * @return
	 */
	public static SparseMatrix solve(SparseMatrixRowMajor A, 
			SparseMatrixRowMajor L, SparseMatrixRowMajor U, SparseMatrix P,
			SparseMatrix X, SparseMatrix F) {
		long begin,end;
		begin = System.currentTimeMillis();
		LU(A, L, U, P);
		end = System.currentTimeMillis();
		System.out.println("FuturEye Sparse LU="+(end-begin)+"ms");
		
		begin = System.currentTimeMillis();
		solvePx(P,X,F);
		end = System.currentTimeMillis();
		System.out.println("FuturEye Sparse SolvePx="+(end-begin)+"ms");
		
		SparseMatrix X2 = new SparseMatrixRowMajor(X.getRowDim(),X.getColDim());
		begin = System.currentTimeMillis();
		solveLx(L,X2,X);
		end = System.currentTimeMillis();
		System.out.println("FuturEye Sparse SolveLx="+(end-begin)+"ms");
		
		begin = System.currentTimeMillis();
		X.clearData();
		solveUx(U,X,X2);
		end = System.currentTimeMillis();
		System.out.println("FuturEye Sparse SolveUx="+(end-begin)+"ms");
		
		return X;
	}
	
	/**
	 * Solve (P*L*U)*X=F with full matrix back substitution,
	 * combined with sparse LU decomposition (A=P*L*U)
	 * this can be faster than any other <tt>solve(...)</tt> function
	 * defined in this class.
	 * 
	 * @param L (Input) Lower triangular matrix with ones on its diagonal from LU decomposition
	 * @param U (Input) Upper triangular matrix from LU decomposition
	 * @param P (Input) Permutation matrix from LU decomposition
	 * @param X (Input) Unknown vectors in matrix form (column based)
	 * @param F (Input) Right hand sides in matrix form (column based)
	 * @return
	 */
	public static FullMatrix backSubstitution(SparseMatrixRowMajor L, SparseMatrixRowMajor U, SparseMatrix P,
			FullMatrix X, FullMatrix F) {
		long begin,end;
		
		begin = System.currentTimeMillis();
		solvePx(P,X,F);
		end = System.currentTimeMillis();
		System.out.println("FuturEye Sparse2 Solve Px="+(end-begin)+"ms");
		
		FullMatrix X2 = X.copy();
		begin = System.currentTimeMillis();
		solveLx(L,X2,X);
		end = System.currentTimeMillis();
		System.out.println("FuturEye Sparse2 Solve Lx="+(end-begin)+"ms");
		
		begin = System.currentTimeMillis();
		solveUx(U,X,X2);
		end = System.currentTimeMillis();
		System.out.println("FuturEye Sparse2 Solve Ux="+(end-begin)+"ms");
		
		return X;
	}
	
	
	public static void test1() {
		SparseMatrix A = new SparseMatrixRowMajor(3,3);
		SparseMatrix L = new SparseMatrixRowMajor(3,3);
		SparseMatrix U = new SparseMatrixRowMajor(3,3);
		SparseMatrix P = new SparseMatrixRowMajor(3,3);
		
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

		SparseVector x = new SparseVectorHashMap(3);
		SparseVector f = new SparseVectorHashMap(3,1.0);
		solve(A,L,U,P,x,f);
		x.print();
		
		SparseVector y = new SparseVectorHashMap(3);
		A.mult(x, y);
		y.print();		
	}
	
	public static void test2() {
		SparseMatrix A = new SparseMatrixRowMajor(3,3);
		SparseMatrix L = new SparseMatrixRowMajor(3,3);
		SparseMatrix U = new SparseMatrixRowMajor(3,3);
		SparseMatrix P = new SparseMatrixRowMajor(3,3);
		
		//double[][] data = {{8,2,9},{4,9,4},{6,7,9}};
		double[][] data = {{1,0,0},{0,0,2},{0,1,-1}};
		for(int i=0;i<data.length;i++) {
			for(int j=0;j<data[i].length;j++)
				A.set(i+1, j+1, data[i][j]);
		}
		A.print();
		
		FullMatrix fA = new FullMatrix(A);
		FullMatrix fL = new FullMatrix(L);
		FullMatrix fU = new FullMatrix(U);
		LU(fA,fL,fU,P);
		fA.print();
		fL.print();
		fU.print();
		P.print();

		FullVector x = new FullVector(3);
		FullVector f = new FullVector(3,1.0);
		solve(fA,fL,fU,P,x,f);
		x.print();
		
		SparseVector y = new SparseVectorHashMap(3);
		A.mult(x.getSparseVector(), y);
		y.print();
	}
	
	public static void test3() {
		SparseMatrix A = new SparseMatrixRowMajor(3,3);
		SparseMatrixRowMajor L = new SparseMatrixRowMajor(3,3);
		SparseMatrixRowMajor U = new SparseMatrixRowMajor(3,3);
		SparseMatrix P = new SparseMatrixRowMajor(3,3);
		
		//double[][] data = {{8,2,9},{4,9,4},{6,7,9}};
		double[][] data = {{1,0,0},{0,0,2},{0,1,-1}};
		for(int i=0;i<data.length;i++) {
			for(int j=0;j<data[i].length;j++)
				A.set(i+1, j+1, data[i][j]);
		}
		A.print();
		
		//LU(A,L,U,P);
		SparseMatrix X = new SparseMatrixRowMajor(3,1);
		SparseMatrix F = new SparseMatrixRowMajor(3,1);
		F.set(1, 1, 1.0);
		F.set(2, 1, 1.0);
		F.set(3, 1, 1.0);
		//solve(A,L,U,P,X,F);
		LU(A,L,U,P);
		FullMatrix XX = new FullMatrix(X);
		FullMatrix FF = new FullMatrix(F);
		backSubstitution(L,U,P,XX,FF);
		A.print();
		L.print();
		U.print();
		P.print();
		XX.print();
		
		
//		SparseMatrix Y = new SparseMatrix(3,1);
//		A.mult(X, Y);
//		Y.print();		

	}
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		test1();
		System.out.println("-----------------------------");
		test2();
		System.out.println("-----------------------------");
		test3();
	}

}
