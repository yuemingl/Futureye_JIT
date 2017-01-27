package edu.uta.futureye.algebra.solver.external;

import org.ejml.alg.dense.decomposition.lu.LUDecompositionAlt_D64;
import org.ejml.alg.dense.linsol.lu.LinearSolverLu_D64;
import org.ejml.data.DenseMatrix64F;

import edu.uta.futureye.algebra.FullMatrix;
import edu.uta.futureye.algebra.SparseMatrixRowMajor;
import edu.uta.futureye.algebra.intf.SparseMatrix;
import edu.uta.futureye.algebra.solver.LUDecomposition;

public class SolverEJML {
	
	/**
	 * Solve <tt>A*X = B</tt> with LU Decomposition.
	 * @param A
	 * @param B
	 * @return X
	 */
	public FullMatrix solve(FullMatrix A,FullMatrix B) {
		DenseMatrix64F AA = new DenseMatrix64F(A.getData());
		DenseMatrix64F BB = new DenseMatrix64F(B.getData());
		
		long begin,end;
		begin = System.currentTimeMillis();
		LUDecompositionAlt_D64 lu = new LUDecompositionAlt_D64();
		lu.decompose(AA);
		end = System.currentTimeMillis();
		System.out.println("EJML LU="+(end-begin)+"ms");
		
		begin = System.currentTimeMillis();
		LinearSolverLu_D64 solver = new LinearSolverLu_D64(lu);
		DenseMatrix64F XX = new DenseMatrix64F(BB.numRows,BB.numCols);
		solver.solve(BB, XX);
		end = System.currentTimeMillis();
		System.out.println("EJML Solve="+(end-begin)+"ms");
		
		FullMatrix rlt = new FullMatrix(BB.numRows,BB.numCols,XX.getData());
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
		SolverEJML solver = new SolverEJML();
		System.out.println("Begin EJML solver...");
		long begin = System.currentTimeMillis();
		FullMatrix X = solver.solve(A, B);
		long end = System.currentTimeMillis();
		System.out.println("EJML time used: "+(end-begin));

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
