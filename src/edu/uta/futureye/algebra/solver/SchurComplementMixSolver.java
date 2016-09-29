package edu.uta.futureye.algebra.solver;

import edu.uta.futureye.algebra.SparseBlockMatrix;
import edu.uta.futureye.algebra.SparseBlockVector;
import edu.uta.futureye.algebra.SparseVectorHashMap;
import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.SparseVector;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.algebra.solver.external.SolverJBLAS;
import edu.uta.futureye.function.FMath;

/**
 * <blockquote><pre>
 * A = (M B')
 *     (B 0 )
 *
 * f = (F)
 *     (G)
 *
 * 求 A*x=f, x=(U P)'
 * 
 * 
 * (B*inv(M)*B')*P = B*inv(M)*F - G
 *             M*U = F - B'*P
 * 
 * where S = B*inv(M)*B' 称为 Schur complement matrix
 *</blockquote></pre>
 *
 * @author liuyueming
 *
 */
public class SchurComplementMixSolver {
	static double eps = 1e-5;
	static double maxIter = 100;

	SparseBlockMatrix A;
	SparseBlockVector f;
	
	public SchurComplementMixSolver(SparseBlockMatrix A, SparseBlockVector f) {
		this.A = A;
		this.f = f;
	}
	
	public SparseBlockVector solve() {
		Matrix M = A.getBlock(1, 1);
		Matrix B = A.getBlock(2, 1);
		Matrix BT = A.getBlock(1, 2);
		Vector F = f.getBlock(1);
		Vector G = f.getBlock(2);
		
		Vector tmp = new SparseVectorHashMap(M.getRowDim());
		Vector rhs = new SparseVectorHashMap(B.getRowDim());
		
		//Schur complement right hand side: B*inv(M)*F - G
		tmp = this.invM_v(F);
		B.mult(tmp, rhs);
		rhs.add(-1.0, G);
		
		// S*P = B*inv(M)*F - G
		// where S = B*inv(M)*B')*P
		Vector P = new SparseVectorHashMap(B.getRowDim(),1.0);
		P = this.CG(rhs, P);
		//P.print();

		//M*U = F - B'*P
		Vector U = new SparseVectorHashMap(BT.getRowDim());
		BT.mult(P, U);
		U.axpy(-1.0, F);
		U = this.invM_v(U);
		
		SparseBlockVector rlt = new SparseBlockVector(2);
		rlt.setBlock(1, (SparseVector)U);
		rlt.setBlock(2, (SparseVector)P);
		return rlt;
		
	}
	
	/**
	 * inv(M)*v
	 * 
	 * @param v
	 * @return
	 */
	public Vector invM_v(Vector v) {
		SolverJBLAS sov = new SolverJBLAS();
		return sov.solveDGESV(A.getBlock(1, 1), v);
	}
	
	/**
	 * S * v = (B*inv(M)*B') * v
	 * @param v
	 * @return
	 */
	public Vector S_v(Vector v) {
		Matrix BT = A.getBlock(1, 2);
		
		Vector tmp = v.copy();
		tmp.setDim(BT.getRowDim());
		BT.mult(v, tmp);
		
		tmp = invM_v(tmp);
		
		Matrix B = A.getBlock(2, 1);
		Vector rlt = v.copy();
		rlt.setDim(BT.getRowDim());
		
		B.mult(tmp, rlt);
		return rlt;
	}
	
	/**
	 * S*x = b
	 * 
	 * S = B*inv(M)*B'
	 * 
	 * @param b
	 * @param x
	 * @return
	 */
    public Vector CG(Vector b, Vector x) {

		double alpha = 0, beta = 0, rho = 0, rho_1 = 0;
		
		Vector r = b.copy(); r.setAll(0.0);
		Vector z = b.copy(); z.setAll(0.0);
		Vector p = b.copy(); p.setAll(0.0);
		Vector q = null;
		
		// r = b - Ax
		//A.multAdd(-1, x, r.set(b));
		r = FMath.axpy(-1.0, this.S_v(x), b);
		
		//for (iter.setFirst(); !iter.converged(r, x); iter.next()) {
		for(int i=0;i<maxIter;i++) {
			double norm2 = r.norm2();
		    System.out.println("Iter----->i="+i+"  norm2="+norm2);
			if(norm2 < eps)
				break;
			
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
		        p.axpy(beta, z);
		    }
		
		    //q = A*p
		    //A.mult(p, q);
		    q= this.S_v(p);
		    alpha = rho / p.dot(q);
		
		    x.add(alpha, p);
		    r.add(-alpha, q);
		
		    rho_1 = rho;
		}
		
		return x;
    }
	
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
