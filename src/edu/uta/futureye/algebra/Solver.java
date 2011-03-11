package edu.uta.futureye.algebra;

import java.util.Map;
import java.util.Map.Entry;

import no.uib.cipr.matrix.sparse.CG;
import no.uib.cipr.matrix.sparse.IterativeSolverNotConvergedException;

import org.netlib.lapack.DGESV;
import org.netlib.util.intW;

import edu.uta.futureye.algebra.intf.AlgebraMatrix;
import edu.uta.futureye.algebra.intf.AlgebraVector;
import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.util.FutureyeException;

public class Solver {
	static double eps = 1e-10;
	
	/**
	 * 
	 * @param m
	 * @param v
	 * @return 方程的解向量，下标从1开始
	 */
	public Vector solve(Matrix m,Vector v) {
		if( !( m.getRowDim() == m.getColDim() &&
				m.getRowDim() == v.getDim()) ) {
			FutureyeException e = new FutureyeException(
					"ERROR: Solver.solver() m.dim!=v.dim ");
			e.printStackTrace();
			System.exit(0);
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
	    
    	Vector rv = new SparseVector(N);
	    for(int i=0;i<N;i++) {
	    	if(Math.abs(b[i][0]) > eps)
	    		rv.set(i+1, b[i][0]);
	    }
	    return rv;

	}
	
	//迭代相对误差
	static double epsIter = 1e-5;
	//迭代绝对误差
	static double epsAbsIter = 1e-15;
	//迭代最大次数
	static double maxIter = 1000;

	public AlgebraVector CG(AlgebraMatrix A, AlgebraVector b, AlgebraVector x) {

		double alpha = 0, beta = 0, rho = 0, rho_1 = 0;
		
		int dim = b.getDim();
		AlgebraVector r = new FullVector(dim);
		AlgebraVector z = new FullVector(dim);
		AlgebraVector p = new FullVector(dim);
		AlgebraVector q = new FullVector(dim);
		
		// r = b - Ax
		//A.multAdd(-1, x, r.set(b));
		A.mult(x, r);
		r.axpy(-1.0, b);
		
		double firstNorm2 = r.norm2();
		double norm2 = 0;
		//for (iter.setFirst(); !iter.converged(r, x); iter.next()) {
		for(int i=0;i<maxIter;i++) {
			norm2 = r.norm2();
			if(norm2<epsIter*firstNorm2) {
				System.out.println("Iter----->i="+i+"  norm2="+norm2);
				return x;
			}
			
		    //M.apply(r, z);
			//Mz=r
			//M：预条件矩阵，取为I,z==r
			z=r;
			
		    rho = r.dot(z);
		
		    if (i==0)
		        p.set(z);
		    else {
		        beta = rho / rho_1;
		        p.axpy(beta, z); //p=beta*p+z
		    }
		
		    //q = A*p
		    A.mult(p, q);
		    alpha = rho / p.dot(q);
		
		    x.add(alpha, p); //x=x+alpha*p
		    r.add(-alpha, q); //r=r-alpha*q
		
		    rho_1 = rho;
		}
		System.out.println("Iter Max----->maxIter="+maxIter+"  norm2="+norm2);
		return x;
    }
	
	public Vector CG1(Matrix A, Vector b, Vector x) {

		double alpha = 0, beta = 0, rho = 0, rho_1 = 0;
		
		
		Vector r = new SparseVector(b.getDim());
		
		Vector z = new SparseVector(b.getDim());
		Vector p = new SparseVector(b.getDim());
		Vector q = new SparseVector(b.getDim());
		
		// r = b - Ax
		//A.multAdd(-1, x, r.set(b));
		A.mult(x, r);
		r = SparseVector.axpy(-1.0, r, b);
		
		double firstNorm2 = r.norm2();
		double norm2 = 0;
		//for (iter.setFirst(); !iter.converged(r, x); iter.next()) {
		for(int i=0;i<maxIter;i++) {
			norm2 = r.norm2();
			if(norm2<epsIter*firstNorm2) {
				System.out.println("----->i="+i+"  norm2="+norm2);
				return x;
			}
			
		    //M.apply(r, z);
			//Mz=r
			//M：预条件矩阵，取为I,z==r
			z=r;
			
		    rho = r.dot(z);
		
		    if (i==0)
		        p.set(z);
		    else {
		        beta = rho / rho_1;
		        //p.scale(beta).add(z);
		        //p=beta*p+z
		        p = SparseVector.axpy(beta, p, z);
		    }
		
		    //q = A*p
		    //A.mult(p, q);
		    A.mult(p, q);
		    alpha = rho / p.dot(q);
		
		    x.add(alpha, p);
		    r.add(-alpha, q);
		    
		    rho_1 = rho;
		}
		System.out.println("Iter Max----->maxIter="+maxIter+"  norm2="+norm2);
		return x;
    }
	
	public Vector CG2(SparseMatrix A, SparseVector b, SparseVector x) {
		int dim = b.getDim();
		no.uib.cipr.matrix.sparse.SparseVector template = 
			new no.uib.cipr.matrix.sparse.SparseVector(dim);
		CG cg = new CG(template);
		
		no.uib.cipr.matrix.sparse.CompRowMatrix A2 = 
			new no.uib.cipr.matrix.sparse.CompRowMatrix(
					A.getRowDim(),A.getColDim(),A.getColIndex());
		
		Map<Integer, Map<Integer, Double>> data = A.getAll();
		for(Entry<Integer,Map<Integer,Double>> row : data.entrySet()) {
			int nRow = row.getKey();
			for(Entry<Integer,Double> col : row.getValue().entrySet()) {
				int nCol = col.getKey();
				A2.set(nRow-1, nCol-1, col.getValue());
			}
			row.getValue().clear();
		}
		data.clear();
		
		no.uib.cipr.matrix.sparse.SparseVector b2 = new no.uib.cipr.matrix.sparse.SparseVector(dim);
		no.uib.cipr.matrix.sparse.SparseVector x2 = new no.uib.cipr.matrix.sparse.SparseVector(dim);
		Map<Integer,Double> bData = b.getAll();
		for(Entry<Integer,Double> ety : bData.entrySet()) {
			b2.set(ety.getKey()-1,ety.getValue());
		}
		for(int i=0;i<dim;i++) {
			x2.set(i, 0.01);
		}
		
		try {
			long begin = System.currentTimeMillis();
			cg.solve(A2, b2, x2);
			long end = System.currentTimeMillis();
			System.out.println("Solve time used:"+(end-begin));
			System.out.println("Iterations: "+cg.getIterationMonitor().iterations());

		} catch (IterativeSolverNotConvergedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		
		for(int i=1;i<=dim;i++) {
			x.set(i, x2.get(i-1));
		}
		
		return x;
    }
	

}
