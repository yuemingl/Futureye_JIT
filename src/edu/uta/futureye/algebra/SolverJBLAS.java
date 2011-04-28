package edu.uta.futureye.algebra;

import org.netlib.lapack.DGESV;
import org.netlib.util.intW;

import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.util.FutureyeException;

public class SolverJBLAS {
	static double eps = 1e-10;
	
	/**
	 * Java BLAS Interface
	 * 
	 * @param m
	 * @param v
	 * @return 方程的解向量，下标从1开始
	 */
	public Vector solveDGESV(Matrix m, Vector v) {
		if( !( m.getRowDim() == m.getColDim() &&
				m.getRowDim() == v.getDim()) ) {
			throw new FutureyeException(
					"ERROR: Solver.solver() m.dim!=v.dim ");
		}
		
	    int N = v.getDim();
	    int nrhs = 1;
	    int[]ipiv = new int[N];
	    
	    double[][]a = new double[N][N];
	    double[][]b = new double[N][1];

	    for(int i=0;i<N;i++) {
	    	for(int j=0;j<N;j++) {
	    		a[i][j] = m.get(i+1, j+1);
	    		//System.out.print(a[i][j]+" ");
	    	}
	    	//System.out.println("");
	    	b[i][0] = v.get(i+1);
	    	//System.out.println(b[i][0]);
	    }
	    intW info = new intW(0);
	    
        System.out.println("Begin Solver...");
        DGESV.DGESV(N, nrhs, a, ipiv, b, info);
        System.out.println("Solver info = " + info.val);
	    
    	Vector rv = v.copy();
	    for(int i=0;i<N;i++) {
	    	if(Math.abs(b[i][0]) > eps)
	    		rv.set(i+1, b[i][0]);
	    }
	    return rv;

	}
	
}
