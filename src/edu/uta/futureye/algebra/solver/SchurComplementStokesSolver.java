/**
 * Copyright (c) 2010, nkliuyueming@gmail.com. All rights reserved.
 * 
 * 
 */
package edu.uta.futureye.algebra.solver;

import edu.uta.futureye.algebra.CompressedColMatrix;
import edu.uta.futureye.algebra.CompressedRowMatrix;
import edu.uta.futureye.algebra.FullMatrix;
import edu.uta.futureye.algebra.FullVector;
import edu.uta.futureye.algebra.SparseBlockMatrix;
import edu.uta.futureye.algebra.SparseBlockVector;
import edu.uta.futureye.algebra.SparseMatrixRowMajor;
import edu.uta.futureye.algebra.SparseVectorHashMap;
import edu.uta.futureye.algebra.intf.AlgebraMatrix;
import edu.uta.futureye.algebra.intf.SparseMatrix;
import edu.uta.futureye.algebra.intf.SparseVector;
import edu.uta.futureye.algebra.solver.external.SolverColt;
import edu.uta.futureye.algebra.solver.external.SolverJBLAS;

/**
 *<blockquote><pre>
 * 2D case:
 * 
 * A = (B1  0   C1)
 *     (0   B2  C2)
 *     (C1' C2' C )
 *
 * x = (u1)
 *     (u2)
 *     (p )
 *
 * f = (f1)
 *     (f2)
 *     (f3)
 *
 * C1' = - trans(C1)
 * C2' = - trans(C2)
 *
 * Solve A*x=f
 * 
 * B1 *u1          + C1*p = f1   ---(1)
 *          B2 *u2 + C2*p = f2   ---(2)
 * C1'*u1 + C2'*u2 +  C*p = f3   ---(3)
 * 
 * From (1)(2), we get
 * u1 = inv(B1) * (f1 - C1*p)  ---(4)
 * u2 = inv(B2) * (f2 - C2*p)  ---(5)
 * 
 * Substitute into (3)：
 * (C - C1'*inv(B1)*C1 - C2'*inv(B2)*C2)*p 
 *                 = f3 - C1'*inv(B1)*f1 - C2'*inv(B2)*f2   ---(6)
 * 
 * where S = C - C1'*inv(B1)*C1 - C2'*inv(B2)*C2 , 称为 Schur complement matrix
 *
 * Solving Steps
 *  Firstly, solve (6), get p
 *  Secondly, from (4)(5), get u1,u2
 *
 *
 * 3D case (@see solve3D()):
 * A = (B1  0   0   C1)
 *     (0   B2  0   C2)
 *     (0   0   B3  C3)
 *     (C1' C2' C3' C )
 *     
 *</blockquote></pre>
 * @author liuyueming
 *
 */
public class SchurComplementStokesSolver {
	protected SparseBlockMatrix A;
	protected SparseBlockVector f;
	double init = 1.0;
	public boolean debug = false;
	
	public SchurComplementStokesSolver(SparseBlockMatrix A,SparseBlockVector f) {
		this.A = A;
		this.f = f;
	}
	
	public void setCGInit(double init) {
		this.init = init;
	}
	
	public SparseBlockVector solve2D() {
		SparseMatrix B1 = A.getBlock(1, 1);
		SparseMatrix B2 = A.getBlock(2, 2);
		SparseMatrix C1 = A.getBlock(1, 3);
		SparseMatrix C2 = A.getBlock(2, 3);
		SparseMatrix C1T = A.getBlock(3, 1);
		SparseMatrix C2T = A.getBlock(3, 2);
		SparseMatrix C = A.getBlock(3, 3);
		SparseVector f1 = f.getBlock(1);
		SparseVector f2 = f.getBlock(2);
		SparseVector f3 = f.getBlock(3);
//		System.out.print("B1=[");B1.print();System.out.print("];");
//		System.out.print("B2=[");B2.print();System.out.print("];");
//		System.out.print("C1=[");C1.print();System.out.print("];");
//		System.out.print("C2=[");C2.print();System.out.print("];");
//		System.out.print("C1T=[");C1T.print();System.out.print("];");
//		System.out.print("C2T=[");C2T.print();System.out.print("];");
//		System.out.print("C=[");C.print();System.out.print("];");
//		System.out.print("f1=[");f1.print();System.out.print("]';");
//		System.out.print("f2=[");f2.print();System.out.print("]';");
//		System.out.print("f3=[");f3.print();System.out.print("]';");
		B1.writeMatFile("B1.mat");
		
		FullVector tmp1 = null;
		FullVector tmp2 = null;
		FullVector rhs  = new FullVector(C1.getColDim());
		FullVector rhs2 = new FullVector(C2.getColDim());
		
//		CompressedRowMatrix BB1 = new CompressedRowMatrix(B1, true);
//		CompressedRowMatrix BB2 = new CompressedRowMatrix(B2, true);
//		CompressedColMatrix CC1 = new CompressedColMatrix(C1, true);
//		CompressedColMatrix CC2 = new CompressedColMatrix(C2, true);
		//for Sparse LU Decomposition
		CompressedRowMatrix BB1 = new CompressedRowMatrix(B1, false);
		CompressedRowMatrix BB2 = new CompressedRowMatrix(B2, false);
		CompressedColMatrix CC1 = new CompressedColMatrix(C1, false);
		CompressedColMatrix CC2 = new CompressedColMatrix(C2, false);
		
		CompressedRowMatrix CC1T = new CompressedRowMatrix(C1T, true);//C1T = - trans(C1)
		CompressedRowMatrix CC2T = new CompressedRowMatrix(C2T, true);//C2T = - trans(C2)
		CompressedRowMatrix CC = new CompressedRowMatrix(C, true);
		FullVector ff1 = new FullVector(f1);
		FullVector ff2 = new FullVector(f2);
		FullVector ff3 = new FullVector(f3);
		
		//Schur complement right hand side: 
		//f3 - C1'*inv(B1)*f1 - C2'*inv(B2)*f2
		tmp1 = invB_v(BB1,ff1);
		tmp2 = invB_v(BB2,ff2);
		CC1T.mult(tmp1, rhs);
		CC2T.mult(tmp2, rhs2);
		rhs.ax(-1.0);
		rhs.add(-1.0, rhs2);
		rhs.add(1.0, ff3);
		
		long begin=0,end=0;
		
//		CompressedColMatrix bc1 = (CompressedColMatrix)invB_C(BB1,CC1);
//		CompressedColMatrix bc2 = (CompressedColMatrix)invB_C(BB2,CC2);
//		CompressedColMatrix bcc1 = (CompressedColMatrix)invB_CColt(BB1,CC1);
//		CompressedColMatrix bcc2 = (CompressedColMatrix)invB_CColt(BB2,CC2);
		
		//S = - C1'*inv(B1)*C1 - C2'*inv(B2)*C2 + C
		CompressedRowMatrix S = new CompressedRowMatrix();
		CompressedRowMatrix S2 = new CompressedRowMatrix();
		
//Dense Colt Solver	
//		CC1T.mult(invB_CColt(BB1,CC1), S);
//		CC2T.mult(invB_CColt(BB2,CC2), S2);
//Sparse Colt Solver		
//		CC1T.mult(invB_CColt(B1,C1), S);
//		CC2T.mult(invB_CColt(B2,C2), S2);		
		
//		CC1T.mult(invB_C(BB1,CC1), S);
//		CC2T.mult(invB_C(BB2,CC2), S2);

//FutuEye LU Decomposition
		CC1T.mult(invB_C(B1,C1), S);
		CC2T.mult(invB_C(B2,C2), S2);
		
//		CC1T.mult(new CompressedColMatrix(invB_C(B1,C1),true), S);
//		CC2T.mult(new CompressedColMatrix(invB_C(B2,C2),true), S2);

		
		S.axpy(-1.0, S2.ax(-1.0)); // S = -S - S2
		S.axpy(1.0, CC);
		
		SolverJBLAS sov = new SolverJBLAS();
		SparseMatrix SS2 = S.getSparseMatrix();
		SparseVector Srhs = rhs.getSparseVector();
		begin = System.currentTimeMillis();
		FullVector p = new FullVector(sov.solveDGESV(SS2, Srhs));
		end = System.currentTimeMillis();
		System.out.println("solveDGESV time="+(end-begin)+"ms");
		
//		Solver sov = new Solver();
//		FullVector p = new FullVector(rhs.getDim(),1.0);
//		sov.solveCG(S, rhs, p);
		
		
		//u1=inv(B1)*(f1-C1*p)
		//u2=inv(B2)*(f2-C2*p)
		FullVector u1 = new FullVector(f1.getDim(),1.0);
		FullVector u2 = new FullVector(f2.getDim(),1.0);
		CC1.getCompressedRowMatrix().mult(p, u1);
		u1.axpy(-1.0, ff1);
		CC2.getCompressedRowMatrix().mult(p, u2);
		u2.axpy(-1.0, ff2);
		u1 = invB_v(BB1,u1);
		u2 = invB_v(BB2,u2);
		
		SparseBlockVector rlt = new SparseBlockVector(3);
		SparseVector uu1 = new SparseVectorHashMap(u1.getData());
		SparseVector uu2 = new SparseVectorHashMap(u2.getData());
		SparseVector pp = new SparseVectorHashMap(p.getData());
		rlt.setBlock(1, uu1);
		rlt.setBlock(2, uu2);
		rlt.setBlock(3, pp);
		return rlt;
	}
	
	/**
	 * <blockquote><pre>
	 * A = (B1  0   0   C1)
	 *     (0   B2  0   C2)
	 *     (0   0   B3  C3)
	 *     (C1' C2' C3' C )
	 *
	 * x = (u1)
	 *     (u2)
	 *     (u3)
	 *     (p )
	 *
	 * f = (f1)
	 *     (f2)
	 *     (f3)
	 *     (f4)
	 *
	 * C1' = - trans(C1)
	 * C2' = - trans(C2)
	 * C3' = - trans(C3)
	 *
	 * Solve A*x=f
	 * 
	 * B1 *u1                   + C1*p = f1   ---(1)
	 *          B2 *u2          + C2*p = f2   ---(2)
	 *                   B3 *u3 + C3*p = f3   ---(3)
	 * C1'*u1 + C2'*u2 + C3'*u3 + C*p  = f3   ---(4)
	 * 
	 * From (1),(2),(3) we get
	 * u1 = inv(B1) * (f1 - C1*p)  ---(5)
	 * u2 = inv(B2) * (f2 - C2*p)  ---(6)
	 * u3 = inv(B3) * (f2 - C2*p)  ---(7)
	 * 
	 * Substitute into (4)：
	 *   (C - C1'*inv(B1)*C1 - C2'*inv(B2)*C2 - C3'*inv(B3)*C3)*p 
	 *  = f4 - C1'*inv(B1)*f1 - C2'*inv(B2)*f2 - C3'*inv(B3)*f3   ---(8)
	 * 
	 * where S = C - C1'*inv(B1)*C1 - C2'*inv(B2)*C2 - C3'*inv(B3)*C3, 称为 Schur complement matrix
	 *
	 * Solving Steps
	 *  Firstly, solve (8), get p
	 *  Secondly, from (5),(6),(7) get u1,u2,u3
	 * </blockquote></pre>
	 * @author liuyueming
	 *
	 */
	public SparseBlockVector solve3D() {
		System.out.println("Begin solve3D...");
		long begin=0,end=0;
		begin = System.currentTimeMillis();
		SparseMatrix B1  = A.getBlock(1, 1);
		SparseMatrix B2  = A.getBlock(2, 2);
		SparseMatrix B3  = A.getBlock(3, 3);
		SparseMatrix C1  = A.getBlock(1, 4);
		SparseMatrix C2  = A.getBlock(2, 4);
		SparseMatrix C3  = A.getBlock(3, 4);
		SparseMatrix C1T = A.getBlock(4, 1);
		SparseMatrix C2T = A.getBlock(4, 2);
		SparseMatrix C3T = A.getBlock(4, 3);
		SparseMatrix C   = A.getBlock(4, 4);
		SparseVector f1 = f.getBlock(1);
		SparseVector f2 = f.getBlock(2);
		SparseVector f3 = f.getBlock(3);
		SparseVector f4 = f.getBlock(4);
		
		FullVector tmp1 = null;
		FullVector tmp2 = null;
		FullVector tmp3 = null;
		FullVector rhs  = new FullVector(C1.getColDim());
		FullVector rhs2 = new FullVector(C2.getColDim());
		FullVector rhs3 = new FullVector(C3.getColDim());
		
//		CompressedRowMatrix BB1  = new CompressedRowMatrix(B1, true);
//		CompressedRowMatrix BB2  = new CompressedRowMatrix(B2, true);
//		CompressedRowMatrix BB3  = new CompressedRowMatrix(B3, true);
//		CompressedColMatrix CC1  = new CompressedColMatrix(C1, true);
//		CompressedColMatrix CC2  = new CompressedColMatrix(C2, true);
//		CompressedColMatrix CC3  = new CompressedColMatrix(C3, true);
		//For sparse LU solver
		CompressedRowMatrix BB1  = new CompressedRowMatrix(B1, false);
		CompressedRowMatrix BB2  = new CompressedRowMatrix(B2, false);
		CompressedRowMatrix BB3  = new CompressedRowMatrix(B3, false);
		CompressedColMatrix CC1  = new CompressedColMatrix(C1, false);
		CompressedColMatrix CC2  = new CompressedColMatrix(C2, false);
		CompressedColMatrix CC3  = new CompressedColMatrix(C3, false);
		
		CompressedRowMatrix CC1T = new CompressedRowMatrix(C1T, true);//C1T = - trans(C1)
		CompressedRowMatrix CC2T = new CompressedRowMatrix(C2T, true);//C2T = - trans(C2)
		CompressedRowMatrix CC3T = new CompressedRowMatrix(C3T, true);//C3T = - trans(C3)
		CompressedRowMatrix CC   = new CompressedRowMatrix(C, true);
		FullVector ff1 = new FullVector(f1);
		FullVector ff2 = new FullVector(f2);
		FullVector ff3 = new FullVector(f3);
		FullVector ff4 = new FullVector(f4);
		end = System.currentTimeMillis();
		System.out.println("Preparing matrix and vector: "+(end-begin)+"ms");
		
		//Schur complement right hand side: 
		//f4 - C1'*inv(B1)*f1 - C2'*inv(B2)*f2 - C3'*inv(B3)*f3
		tmp1 = invB_v(BB1,ff1);
		tmp2 = invB_v(BB2,ff2);
		tmp3 = invB_v(BB3,ff3);
		CC1T.mult(tmp1, rhs);
		CC2T.mult(tmp2, rhs2);
		CC3T.mult(tmp3, rhs3);
		rhs.ax(-1.0);
		rhs.add(-1.0, rhs2);
		rhs.add(-1.0, rhs3);
		rhs.add(1.0, ff4);
		
		

		
		
		//S = C - C1'*inv(B1)*C1 - C2'*inv(B2)*C2 - C3'*inv(B3)*C3
		CompressedRowMatrix S = new CompressedRowMatrix();
		CompressedRowMatrix S2 = new CompressedRowMatrix();
		CompressedRowMatrix S3 = new CompressedRowMatrix();
		
//		CC1T.mult(invB_CColt(BB1,CC1), S);
//		CC2T.mult(invB_CColt(BB2,CC2), S2);
//		CC3T.mult(invB_CColt(BB3,CC3), S3);
		
//		CC1T.mult(invB_C(BB1,CC1), S);
//		CC2T.mult(invB_C(BB2,CC2), S2);
//		CC3T.mult(invB_C(BB3,CC3), S3);
		
		CC1T.mult(invB_C(B1,C1), S);
		CC2T.mult(invB_C(B2,C2), S2);
		CC3T.mult(invB_C(B3,C3), S3);
		
		S.axpy(-1.0, S2.ax(-1.0));//S  = -S -S2
		S.axpy(1.0, S3.ax(-1.0)); //S += -S3
		S.axpy(1.0, CC);
		
		
		SolverJBLAS sov = new SolverJBLAS();
		SparseMatrix SS2 = S.getSparseMatrix();
		SparseVector Srhs = rhs.getSparseVector();
		begin = System.currentTimeMillis();
		FullVector p = new FullVector(sov.solveDGESV(SS2, Srhs));
		end = System.currentTimeMillis();
		System.out.println("solveDGESV time="+(end-begin)+"ms");
		
//		Solver sov = new Solver();
//		FullVector p = new FullVector(rhs.getDim(),1.0);
//		sov.solveCG(S, rhs, p);
		
		
		//u1=inv(B1)*(f1-C1*p)
		//u2=inv(B2)*(f2-C2*p)
		FullVector u1 = new FullVector(f1.getDim(),1.0);
		FullVector u2 = new FullVector(f2.getDim(),1.0);
		FullVector u3 = new FullVector(f3.getDim(),1.0);
		CC1.getCompressedRowMatrix().mult(p, u1);
		u1.axpy(-1.0, ff1);
		CC2.getCompressedRowMatrix().mult(p, u2);
		u2.axpy(-1.0, ff2);
		CC3.getCompressedRowMatrix().mult(p, u3);
		u3.axpy(-1.0, ff3);
		u1 = invB_v(BB1,u1);
		u2 = invB_v(BB2,u2);
		u3 = invB_v(BB3,u3);
		
		SparseBlockVector rlt = new SparseBlockVector(4);
		SparseVector uu1 = new SparseVectorHashMap(u1.getData());
		SparseVector uu2 = new SparseVectorHashMap(u2.getData());
		SparseVector uu3 = new SparseVectorHashMap(u3.getData());
		SparseVector pp  = new SparseVectorHashMap(p.getData());
		rlt.setBlock(1, uu1);
		rlt.setBlock(2, uu2);
		rlt.setBlock(3, uu3);
		rlt.setBlock(4, pp);
		return rlt;
	}
	
	
	/**
	 * inv(B)*v
	 * 
	 * @param v
	 * @return
	 */
	public FullVector invB_v(CompressedRowMatrix B, FullVector v) {
//		FullVector x = new FullVector(v.getDim(),0.1);
//		x.setRandom(init,0.0);
//		Solver sov = new Solver();
//		sov.debug = this.debug;
//		sov.solveCGS(B, v, x);
		
		SolverJBLAS sov = new SolverJBLAS();
		SparseMatrix SS2 = B.getSparseMatrix();
		SparseVector Srhs = v.getSparseVector();
		long begin=0,end=0;
		begin = System.currentTimeMillis();
		FullVector x = new FullVector(sov.solveDGESV(SS2, Srhs));
		end = System.currentTimeMillis();
		System.out.println("solve inv(B)*v time="+(end-begin)+"ms");
		return x;
	}
	
	public FullMatrix invB_C(SparseMatrix B, SparseMatrix C) {
		long begin,end;
		begin = System.currentTimeMillis();
		SparseMatrixRowMajor L = new SparseMatrixRowMajor(B.getRowDim(),B.getColDim());
		SparseMatrixRowMajor U = new SparseMatrixRowMajor(B.getRowDim(),B.getColDim());
		SparseMatrix P = new SparseMatrixRowMajor(B.getRowDim(),B.getColDim());
		SparseMatrix X = new SparseMatrixRowMajor(C.getRowDim(),C.getColDim());
//		LUDecomposition.solve((SparseMatrixRowMajor)B, L, U, P, X, C);
		LUDecomposition.LU((SparseMatrixRowMajor)B, L, U, P);
//		FullMatrix LL = new FullMatrix(L);
//		FullMatrix UU = new FullMatrix(U);
		FullMatrix CC = new FullMatrix(C);
		FullMatrix XX = new FullMatrix(X);
		LUDecomposition.backSubstitution(L, U, P, XX, CC);
		end = System.currentTimeMillis();
		System.out.println("SparseLU Solve time="+(end-begin)+"ms");
		return XX;
	}
	
	/**
	 * inv(B)*C
	 * @param B
	 * @param C
	 * @return
	 */
	public AlgebraMatrix invB_C(CompressedRowMatrix B, CompressedColMatrix C) {
		CompressedColMatrix BC = new CompressedColMatrix(B.getRowDim(),C.getColDim());
		
		FullVector v = new FullVector(C.getRowDim());
		FullVector x = new FullVector(C.getRowDim());
		
//		x.setRandom(init,0.0);
//		Solver sov = new Solver();
//		sov.debug = this.debug;
//		int colDim = C.getColDim();
//		for(int c=1; c<=colDim; c++) {
//			C.getColVector(c, v);
//			sov.solveCGS(B, v, x);
//			FullVector.SparseData sd = x.getSparseData();
//			BC.setCol(c, sd.index, sd.data);
//		}
		
//		SolverMTJ solMTJ = new SolverMTJ();
//		solMTJ.debug = this.debug;
//		int colDim = C.getColDim();
//		for(int c=1; c<=colDim; c++) {
//			C.getColVector(c, v);
//			SparseVector xx = x.getSparseVector();
//			//solMTJ.solveGMRES(B.getSparseMatrix(), v.getSparseVector(), xx);
//			//solMTJ.solveBiCG(B.getSparseMatrix(), v.getSparseVector(), xx);
//			//solMTJ.solveBiCGstab(B.getSparseMatrix(), v.getSparseVector(), xx);
//			//solMTJ.solveQMR(B.getSparseMatrix(), v.getSparseVector(), xx);
//			solMTJ.solveCG(B.getSparseMatrix(), v.getSparseVector(), xx);
//			//solMTJ.solveCGS(B.getSparseMatrix(), v.getSparseVector(), xx);
//			FullVector.SparseData sd = new FullVector(xx).getSparseData();
//			BC.setCol(c, sd.index, sd.data);
//		}
		
		
//		//LU分解SparseMatrix
//		int N = B.getRowDim();
//		SparseMatrix A = B.getSparseMatrix();
//		SparseMatrix L = new SparseMatrix(N,N);
//		SparseMatrix U = new SparseMatrix(N,N);
//		SparseMatrix P = new SparseMatrix(N,N);
//		
//		SparseVector x = new SparseVector(N);
//		Vector x2 = x.copy();
//		LUDecomposition.LU(A, L, U, P);
//		
//		int colDim = C.getColDim();
//		FullVector v = new FullVector(C.getRowDim());
//		for(int c=1; c<=colDim; c++) {
//			System.out.println(c+"/"+colDim);
//			C.getColVector(c, v);
//			Vector f = v.getSparseVector();
//			LUDecomposition.solvePx(P,x,f);
//			LUDecomposition.solveLx(L,x2,x);
//			LUDecomposition.solveUx(U,x,x2);
//			FullVector fx = new FullVector(x);
//			FullVector.SparseData sd = fx.getSparseData();
//			BC.setCol(c, sd.index, sd.data);
//		}
		
		//LU分解FullMatrix
		int N = B.getRowDim();
		FullMatrix fA = new FullMatrix(B.getSparseMatrix());
		FullMatrix fL = new FullMatrix(N,N);
		FullMatrix fU = new FullMatrix(N,N);
		SparseMatrix P = new SparseMatrixRowMajor(N,N);
		
		long begin=0,end=0;
//		FullVector x2 = x.copy();
//		begin = System.currentTimeMillis();
//		LUDecomposition.LU(fA, fL, fU, P);
//		end = System.currentTimeMillis();
//		System.out.println("LUDecomposition time="+(end-begin)+"ms");
//		
//		int colDim = C.getColDim();
//		begin = System.currentTimeMillis();
//		for(int c=1; c<=colDim; c++) {
//			//System.out.println(c+"/"+colDim);
//			C.getColVector(c, v);
//			LUDecomposition.solvePx(P,x,v);
//			LUDecomposition.solveLx(fL,x2,x);
//			LUDecomposition.solveUx(fU,x,x2);
//			FullVector.SparseData sd = x.getSparseData();
//			BC.setCol(c, sd.index, sd.data);
//		}
//		end = System.currentTimeMillis();
//		System.out.println("Back substitution time="+(end-begin)+"ms");
//		return BC;
		begin = System.currentTimeMillis();
		FullMatrix fC = new FullMatrix(C.getSparseMatrix());
		FullMatrix fX = new FullMatrix(C.getRowDim(),C.getColDim());
		LUDecomposition.solve(fA, fL, fU, P, fX, fC);
		end = System.currentTimeMillis();
		System.out.println("LUDecomposition.solve() time="+(end-begin)+"ms");
		return fX;

		
//		FullMatrix fB = new FullMatrix(B.getSparseMatrix());
//		FullMatrix fC = new FullMatrix(C.getSparseMatrix());
//		SolverColt solver = new SolverColt();
//		long begin=0,end=0;
//		begin = System.currentTimeMillis();
//		FullMatrix fBC = solver.solve(fB, fC);
//		int[] colIdx = new int[fBC.getRowDim()];
//		double[] col = new double[fBC.getRowDim()];
//		for(int c=0;c<fBC.getColDim();c++) {
//			for(int r=0;r<fBC.getRowDim();r++) {
//				colIdx[r] = r;
//				col[r] = fBC.data[r][c];
//			}
//			BC.setCol(c+1, colIdx, col);
//		}
//		end = System.currentTimeMillis();
//		System.out.println("SolverColt time="+(end-begin)+"ms");
//		return BC;
	}
	
	public AlgebraMatrix invB_CColt(CompressedRowMatrix B, CompressedColMatrix C) {
//		CompressedColMatrix BC = new CompressedColMatrix(B.rowDim,C.colDim);
		
		FullMatrix fB = new FullMatrix(B.getSparseMatrix());
		FullMatrix fC = new FullMatrix(C.getSparseMatrix());
		SolverColt solver = new SolverColt();
		long begin=0,end=0;
		begin = System.currentTimeMillis();
		FullMatrix fBC = solver.solve(fB, fC);
//		int[] colIdx = new int[fBC.getRowDim()];
//		double[] col = new double[fBC.getRowDim()];
//		for(int c=0;c<fBC.getColDim();c++) {
//			for(int r=0;r<fBC.getRowDim();r++) {
//				colIdx[r] = r;
//				col[r] = fBC.data[r][c];
//			}
//			BC.setCol(c+1, colIdx, col);
//		}
		end = System.currentTimeMillis();
		System.out.println("Dense Solver: Colt time="+(end-begin)+"ms");
		
		return fBC;
	}

	public AlgebraMatrix invB_CColt(SparseMatrix B, SparseMatrix C) {
		SolverColt solver = new SolverColt();
		long begin=0,end=0;
		begin = System.currentTimeMillis();
		FullMatrix fBC = solver.solve(B, C);
		end = System.currentTimeMillis();
		System.out.println("Sparse Solver: Colt time="+(end-begin)+"ms");
		
		return fBC;
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}
}
