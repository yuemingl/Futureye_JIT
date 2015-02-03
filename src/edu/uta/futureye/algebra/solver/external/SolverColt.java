package edu.uta.futureye.algebra.solver.external;

import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.impl.SparseDoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import edu.uta.futureye.algebra.FullMatrix;
import edu.uta.futureye.algebra.FullVector;
import edu.uta.futureye.algebra.SparseMatrixRowMajor;
import edu.uta.futureye.algebra.intf.MatrixEntry;
import edu.uta.futureye.algebra.intf.SparseMatrix;
import edu.uta.futureye.algebra.solver.LUDecomposition;

/**
 * Colt Interface
 * 
 * @author liuyueming
 *
 */
public class SolverColt {
	
	/**
	 * Solve <tt>A*X = B</tt> with LU Decomposition.
	 * @param A
	 * @param B
	 * @return X
	 */
	public FullMatrix solve(FullMatrix A,FullMatrix B) {
		
		DoubleMatrix2D AA = new DenseDoubleMatrix2D(A.getData());
		DoubleMatrix2D BB = new DenseDoubleMatrix2D(B.getData());
//		Algebra a = new Algebra();
//		DoubleMatrix2D X = a.solve(AA, BB);
		long begin,end;
		begin = System.currentTimeMillis();
		cern.colt.matrix.linalg.LUDecomposition lu = 
			new cern.colt.matrix.linalg.LUDecomposition(AA);
		end = System.currentTimeMillis();
		System.out.println("Colt LU="+(end-begin)+"ms");
		
		begin = System.currentTimeMillis();
		DoubleMatrix2D X = lu.solve(BB);
		end = System.currentTimeMillis();
		System.out.println("Colt Solve="+(end-begin)+"ms");
		
		FullMatrix rlt = new FullMatrix(X.toArray(),true);
		return rlt;
	}
	
	/**
	 * Solve <tt>A*X = B</tt> with LU Decomposition.
	 * @param A
	 * @param B
	 * @return X
	 */
	public FullMatrix solve(SparseMatrix A,SparseMatrix B) {
		DoubleMatrix2D AA = new SparseDoubleMatrix2D(A.getRowDim(),A.getColDim());
		DoubleMatrix2D BB = new SparseDoubleMatrix2D(B.getRowDim(),B.getColDim());
		for(MatrixEntry e : A) 
			AA.setQuick(e.getRow()-1, e.getCol()-1, e.getValue());
		for(MatrixEntry e : B) 
			BB.setQuick(e.getRow()-1, e.getCol()-1, e.getValue());
		
//		Algebra a = new Algebra();
//		DoubleMatrix2D X = a.solve(AA, BB);
		long begin,end;
		begin = System.currentTimeMillis();
		cern.colt.matrix.linalg.LUDecomposition lu = 
			new cern.colt.matrix.linalg.LUDecomposition(AA);
		end = System.currentTimeMillis();
		System.out.println("Sparse Colt LU="+(end-begin)+"ms");
		
		begin = System.currentTimeMillis();
		DoubleMatrix2D X = lu.solve(BB);
		end = System.currentTimeMillis();
		System.out.println("Sparse Colt Solve="+(end-begin)+"ms");
		
		FullMatrix rlt = new FullMatrix(X.toArray(),true);
		return rlt;
	}
	
	
	/**
	 * Solve <tt>A*x = b</tt> with LU Decomposition.
	 * @param A
	 * @param b
	 * @return
	 */
	public FullVector solve(FullMatrix A,FullVector b) {
		
		DoubleMatrix2D mat = new DenseDoubleMatrix2D(A.getData());
		cern.colt.matrix.linalg.LUDecomposition lu = 
			new cern.colt.matrix.linalg.LUDecomposition(mat);
		double [][]dBM=new double[b.getDim()][1];
		double[] dB = b.getData();
		for(int i=0;i<b.getDim();i++) {
			dBM[i][0] = dB[i];
		}
		DoubleMatrix2D X = lu.solve(new DenseDoubleMatrix2D(dBM));
		
		double[] dR = new double[b.getDim()];
		for(int i=0;i<b.getDim();i++) {
			dR[i] = X.get(i, 0);
		}
		FullVector rlt = new FullVector(dR,false);
		return rlt;
	}
	
	public static void test4() {
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
		SolverColt solver = new SolverColt();
		System.out.println("Begin Colt solver...");
		long begin = System.currentTimeMillis();
		FullMatrix X = solver.solve(A, B);
		long end = System.currentTimeMillis();
		System.out.println("Colt time used: "+(end-begin));

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
	
	public static void test3() {
		double[][] dataA = {{8,2,9},{4,9,4},{6,7,9}};
		double[][] dataB = {{1,0,0},{0,0,2},{0,1,-1}};
		FullMatrix A = new FullMatrix(dataA,true);
		FullMatrix B = new FullMatrix(dataB,true);
		SolverColt solver = new SolverColt();
		FullMatrix X = solver.solve(A, B);
		A.print();
		B.print();
		X.print();
		A.mult(X, B);
		B.print();
	}
	
	public static void test2() {
		double[][] data = {{8,2,9},{4,9,4},{6,7,9}};
		//double[][] data = {{1,0,0},{0,0,2},{0,1,-1}};

		DoubleMatrix2D mat = new DenseDoubleMatrix2D(data);
		cern.colt.matrix.linalg.LUDecomposition lu = new cern.colt.matrix.linalg.LUDecomposition(mat);
		System.out.println(lu.getL().toString());
	}
	
	public static void test1() {
       // For Colt
       DoubleMatrix2D xmatrix,ymatrix,zmatrix;

       DoubleFactory2D myfactory;
       myfactory = DoubleFactory2D.dense;
       xmatrix = myfactory.make(8,2);
       ymatrix = myfactory.make(8,1);

       xmatrix.set(0,0,1);
       xmatrix.set(1,0,1);
       xmatrix.set(2,0,1);
       xmatrix.set(3,0,1);
       xmatrix.set(4,0,1);
       xmatrix.set(5,0,1); 
       xmatrix.set(6,0,1);
       xmatrix.set(7,0,1);

       xmatrix.set(0,1,80);
       xmatrix.set(1,1,220);
       xmatrix.set(2,1,140);
       xmatrix.set(3,1,120);
       xmatrix.set(4,1,180);
       xmatrix.set(5,1,100);
       xmatrix.set(6,1,200);
       xmatrix.set(7,1,160);

       ymatrix.set(0,0,0.6);
       ymatrix.set(1,0,6.70);
       ymatrix.set(2,0,5.30);
       ymatrix.set(3,0,4.00);
       ymatrix.set(4,0,6.55);
       ymatrix.set(5,0,2.15);
       ymatrix.set(6,0,6.60);
       ymatrix.set(7,0,5.75);

       Algebra myAlgebra = new Algebra();
       zmatrix = myAlgebra.solve(xmatrix,ymatrix);
       System.err.println(xmatrix);
       System.err.println(ymatrix);
       System.err.println(zmatrix);	
	}
	
	
	public static void main(String[] args) {
		//test1();
		//test2();
		//test3();
		test4();
	}
	
}
