package edu.uta.futureye.algebra;

import edu.uta.futureye.algebra.intf.BlockMatrix;
import edu.uta.futureye.algebra.intf.BlockVector;

/**
 * A = (B1  0   C1)
 *     (0   B2  C2)
 *     (C1' C2' C )
 *
 * f = (f1)
 *     (f2)
 *     (f3)
 *     
 * x = (u1)
 *     (u2)
 *     (p )
 *     
 * C1' = - trans(C1)
 * C2' = - trans(C2)
 *
 * 求 A*x=f, x=(u1 u2 p)'
 * 
 * B1 *u1 + C1 *p  = f1   ---(1)
 * B2 *u2 + C2 *p  = f2   ---(2)
 * C1'*u1 + C2'*u2 = f3   ---(3)
 * 
 * 由(1)(2)得：
 * u1 = inv(B1) * (f1 - C1*p)  ---(4)
 * u2 = inv(B2) * (f2 - C2*p)  ---(5)
 * 
 * 带入(3)：
 * (C1'*inv(B1)*C1 + C2'*inv(B2)*C2 + C)*p 
 *                 = C1'*inv(B1)*f1 + C2'*inv(B2)*f2 + f3   ---(6)
 *                 
 * where S = C1'*inv(B1)*C1 + C2'*inv(B2)*C2 + C, 称为 Schur complement matrix
 *
 *求解：
 *  先求解(6)得p，再由(4)(5)得u1,u2
 *
 * @author liuyueming
 *
 */
public class SchurComplementStokesSolver {
	protected BlockMatrix A;
	protected BlockVector f;
	
	public SchurComplementStokesSolver(BlockMatrix A,BlockVector f) {
		this.A = A;
		this.f = f;
	}
	
	public BlockVector solve() {
		SparseMatrix B1 = (SparseMatrix)A.getBlock(1, 1);
		SparseMatrix B2 = (SparseMatrix)A.getBlock(2, 2);
		SparseMatrix C1 = (SparseMatrix)A.getBlock(1, 3);
		SparseMatrix C2 = (SparseMatrix)A.getBlock(2, 3);
		SparseMatrix C1T = (SparseMatrix)A.getBlock(3, 1);
		SparseMatrix C2T = (SparseMatrix)A.getBlock(3, 2);
		SparseMatrix C = (SparseMatrix)A.getBlock(3, 3);
		SparseVector f1 = (SparseVector)f.getBlock(1);
		SparseVector f2 = (SparseVector)f.getBlock(2);
		SparseVector f3 = (SparseVector)f.getBlock(3);
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
		
		FullVector tmp1 = null;
		FullVector tmp2 = null;
		FullVector rhs  = new FullVector(C1.getColDim());
		FullVector rhs2 = new FullVector(C2.getColDim());
		
		CompressedRowMatrix BB1 = new CompressedRowMatrix(B1, true);
		CompressedRowMatrix BB2 = new CompressedRowMatrix(B2, true);
		CompressedColMatrix CC1 = new CompressedColMatrix(C1, true);
		CompressedColMatrix CC2 = new CompressedColMatrix(C2, true);
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
		
		//S = - C1'*inv(B1)*C1 - C2'*inv(B2)*C2 + C
		CompressedRowMatrix S = new CompressedRowMatrix();
		CompressedRowMatrix S2 = new CompressedRowMatrix();
		CC1T.mult(invB_C(BB1,CC1), S);
		CC2T.mult(invB_C(BB2,CC2), S2);
		S.axpy(-1.0, S2.ax(-1.0));
		S.axpy(1.0, CC);
		
		Solver sov = new Solver();
		FullVector p = new FullVector(rhs.getDim(),1.0);
		sov.CG(S, rhs, p);
		
		
		//u1=inv(B1)*(f1-C1*p)
		//u2=inv(B2)*(f2-C2*p)
		FullVector u1 = new FullVector(f1.getDim(),1.0);
		FullVector u2 = new FullVector(f2.getDim(),1.0);
		CC1.convertToCompressedRow().mult(p, u1);
		u1.axpy(-1.0, ff1);
		CC2.convertToCompressedRow().mult(p, u2);
		u2.axpy(-1.0, ff2);
		u1 = invB_v(BB1,u1);
		u2 = invB_v(BB2,u2);
		
		SparseBlockVector rlt = new SparseBlockVector(3);
		SparseVector uu1 = new SparseVector(u1.getData());
		SparseVector uu2 = new SparseVector(u2.getData());
		SparseVector pp = new SparseVector(p.getData());
		rlt.setBlock(1, uu1);
		rlt.setBlock(2, uu2);
		rlt.setBlock(3, pp);
		return rlt;
	}
	
	/**
	 * inv(B)*v
	 * 
	 * @param v
	 * @return
	 */
	public FullVector invB_v(CompressedRowMatrix B, FullVector v) {
		FullVector x = new FullVector(v.getDim(),0.1);
		x.setRandom(1.0);
		Solver sov = new Solver();
		sov.CGS(B, v, x);
		return x;
	}
	
	/**
	 * inv(B)*C
	 * @param B
	 * @param C
	 * @return
	 */
	public CompressedColMatrix invB_C(CompressedRowMatrix B, CompressedColMatrix C) {
		CompressedColMatrix BC = new CompressedColMatrix(B.rowDim,C.colDim);
		FullVector v = new FullVector(C.getRowDim());
		FullVector x = new FullVector(C.getRowDim(),1.0);
		
		Solver sov = new Solver();
		int colDim = C.getColDim();
		for(int c=1; c<=colDim; c++) {
			C.getColVector(c, v);
			sov.CGS(B, v, x);
			FullVector.SparseData sd = x.getSparseData();
			BC.setCol(c, sd.index, sd.data);
		}
		return BC;
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}
}
