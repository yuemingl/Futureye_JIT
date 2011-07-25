package edu.uta.futureye.algebra;

import edu.uta.futureye.algebra.intf.AlgebraMatrix;
import edu.uta.futureye.algebra.intf.AlgebraVector;
import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.util.FutureyeException;

public class Solver {
	//迭代相对误差
	public double epsIter = 1e-5;
	//迭代绝对误差
	public double epsAbsIter = 1e-15;
	//迭代最大次数
	public double maxIter = 10000;
	
	/**
	 * Conjugate Gradients iterative method, solves 
	 * symmetric positive definite linear system:
	 * <tt>Ax = b</tt>
	 * 
	 * @param A
	 * @param b
	 * @param x
	 * @return
	 */
	public AlgebraVector solveCG(AlgebraMatrix A, AlgebraVector b, 
			AlgebraVector x) {

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
			if(norm2<=this.epsIter*firstNorm2 || norm2<=this.epsAbsIter) {
				System.out.println(
						String.format("Iter----->i=%05d, RError=%8.3e, AError=%8.3e", 
								i,norm2/firstNorm2,norm2));
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
	
	/**
	 * Conjugate Gradients squared iterative method,
	 * solves the unsymmetric linear system
	 * <tt>Ax = b</tt>
	 * 
	 * @param A
	 * @param b
	 * @param x
	 * @return
	 */
	public AlgebraVector solveCGS(AlgebraMatrix A, AlgebraVector b, 
			AlgebraVector x) {

        double rho_1 = 0, rho_2 = 0, alpha = 0, beta = 0;
		
		int dim = b.getDim();
		AlgebraVector r = new FullVector(dim);
		AlgebraVector p = new FullVector(dim);
		AlgebraVector q = new FullVector(dim);
		AlgebraVector u = new FullVector(dim);
		AlgebraVector phat = new FullVector(dim);
		AlgebraVector qhat = new FullVector(dim);
		AlgebraVector vhat = new FullVector(dim);
		AlgebraVector uhat = new FullVector(dim);
		AlgebraVector sum = new FullVector(dim);
		AlgebraVector rtilde = new FullVector(dim);
        
		// r = b - Ax
		//A.multAdd(-1, x, r.set(b));
		A.mult(x, r);
		r.axpy(-1.0, b);
	    rtilde.set(r);
		
		double firstNorm2 = r.norm2();
		double norm2 = 0;
		//for (iter.setFirst(); !iter.converged(r, x); iter.next()) {
		for(int i=0;i<maxIter;i++) {
			norm2 = r.norm2();
			if(norm2<=this.epsIter*firstNorm2 || norm2<=this.epsAbsIter) {
				System.out.println(
						String.format("Iter----->i=%05d, RError=%8.3e, AError=%8.3e", 
								i,norm2/firstNorm2,norm2));
				return x;
			}
			
            rho_1 = rtilde.dot(r);
            if (rho_1 == 0) {
        		//System.out.println("Iter NotConverge maxIter="+i+"  norm2="+norm2);
            	//return x;
        		throw new FutureyeException("NotConverge, rho_1==0, iter="+i);
            }

            if (i==0) {
                u.set(r);
                p.set(u);
            } else {
                beta = rho_1 / rho_2;
                u.set(r).add(beta, q);
                sum.set(q).add(beta, p);
                p.set(u).add(beta, sum);
            }

            //M.apply(p, phat);
            phat.set(p);
            
            A.mult(phat, vhat);
            alpha = rho_1 / rtilde.dot(vhat);
            q.set(-alpha, vhat).add(u);
            
            //M.apply(sum.set(u).add(q), uhat);
            uhat.set(sum.set(u).add(q));
            x.add(alpha, uhat);
            A.mult(uhat, qhat);
            r.add(-alpha, qhat);

            rho_2 = rho_1;
        }
		System.out.println("Iter Max----->maxIter="+maxIter+"  norm2="+norm2);
		return x;
    }	
	
	/////////////////////////////////////////////////////////////
	
	public Vector solveCG(Matrix A, Vector b, Vector x) {
		if( !( A.getRowDim() == A.getColDim() &&
				A.getRowDim() == b.getDim()) ) {
			throw new FutureyeException(
					"ERROR: Solver.solver() m.dim!=v.dim ");
		}
		AlgebraMatrix algStiff = new CompressedRowMatrix((SparseMatrix)A,false);
		FullVector algLoad = new FullVector(b);
		FullVector algU = new FullVector(x);
		solveCG(algStiff, algLoad, algU);
		double[] data = algU.getData();
		for(int i=0;i<data.length;i++) {
			x.set(i+1, data[i]);
		}
		return x;
	}
	
	public Vector solveCG(Matrix A, Vector b) {
		SparseVector x = new SparseVector(b.getDim(),0.1);
		return solveCG(A,b,x);
	}
	
	public Vector solveCGS(Matrix A, Vector b, Vector x) {
		if( !( A.getRowDim() == A.getColDim() &&
				A.getRowDim() == b.getDim()) ) {
			throw new FutureyeException(
					"ERROR: Solver.solver() m.dim!=v.dim ");
		}
		//CGS
		AlgebraMatrix algStiff = new CompressedRowMatrix((SparseMatrix)A,false);
		FullVector algLoad = new FullVector(b);
		FullVector algU = new FullVector(x);
		solveCGS(algStiff, algLoad, algU);
		double[] data = algU.getData();
		for(int i=0;i<data.length;i++) {
			x.set(i+1, data[i]);
		}
		return x;
	}
	
	public Vector solveCGS(Matrix A, Vector b) {
		Vector x = b.copy();
		for(int i=1;i<x.getDim();i++) {
			x.set(i, 0.1);
		}
		return solveCGS(A,b,x);
	}
	
	////////////////////////////////////////////////////////////////
	
	@Deprecated
	public Vector solveCG2(Matrix A, Vector b, Vector x) {

		double alpha = 0, beta = 0, rho = 0, rho_1 = 0;
		
		
		Vector r = new SparseVector(b.getDim());
		
		Vector z = new SparseVector(b.getDim());
		Vector p = new SparseVector(b.getDim());
		Vector q = new SparseVector(b.getDim());
		
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
		        p.axpy(beta, z);
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
	

}
