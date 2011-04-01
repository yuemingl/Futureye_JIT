package edu.uta.futureye.algebra;

import edu.uta.futureye.algebra.intf.BlockMatrix;
import edu.uta.futureye.algebra.intf.BlockVector;
import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.Vector;

/**
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
 *
 * @author liuyueming
 *
 */
public class SchurComplementSolver {
	static double eps = 1e-5;
	static double maxIter = 100;

	BlockMatrix A;
	BlockVector f;
	
	public SchurComplementSolver(BlockMatrix A,BlockVector f) {
		this.A = A;
		this.f = f;
	}
	
	public Vector solve() {
		Matrix M = A.getBlock(1, 1);
		Matrix B = A.getBlock(2, 1);
		Matrix BT = A.getBlock(1, 2);
		Vector F = f.getBlock(1);
		Vector G = f.getBlock(2);
		
		Vector tmp = new SparseVector(M.getRowDim());
		Vector rhs = new SparseVector(B.getRowDim());
		
		//Schur complement right hand side: B*inv(M)*F - G
		tmp = this.invM_v(F);
		B.mult(tmp, rhs);
		rhs = SparseVector.axpy(-1.0, G, rhs);
		
		// S*P = B*inv(M)*F - G
		// where S = B*inv(M)*B')*P
		Vector P = new SparseVector(B.getRowDim(),1.0);
		P = this.CG(rhs, P);
		//P.print();

		//M*U = F - B'*P
		Vector U = new SparseVector(BT.getRowDim());
		BT.mult(P, U);
		U = SparseVector.axpy(-1.0, U, F);
		U = this.invM_v(U);
		
		SparseBlockVector rlt = new SparseBlockVector(2);
		rlt.setBlock(1, U);
		rlt.setBlock(2, P);
		return rlt;
		
	}
	
	/**
	 * inv(M)*v
	 * 
	 * @param v
	 * @return
	 */
	public Vector invM_v(Vector v) {
		Solver sov = new Solver();
		return sov.solve(A.getBlock(1, 1), v);
	}
	
	/**
	 * S * v = (B*inv(M)*B') * v
	 * @param v
	 * @return
	 */
	public Vector S_v(Vector v) {
		Matrix BT = A.getBlock(1, 2);
		Vector tmp = new SparseVector(BT.getRowDim());
		BT.mult(v, tmp);
		
		tmp = invM_v(tmp);
		
		Matrix B = A.getBlock(2, 1);
		Vector rlt = new SparseVector(B.getRowDim());
		
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
		
		Vector r = new SparseVector(b.getDim());
		
		Vector z = new SparseVector(b.getDim());
		Vector p = new SparseVector(b.getDim());
		Vector q = null;
		
		// r = b - Ax
		//A.multAdd(-1, x, r.set(b));
		r = SparseVector.axpy(-1.0, this.S_v(x), b);
		
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
		        p = SparseVector.axpy(beta, p, z);
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
