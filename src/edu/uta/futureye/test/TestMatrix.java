package edu.uta.futureye.test;

import edu.uta.futureye.algebra.CompressedColMatrix;
import edu.uta.futureye.algebra.CompressedRowMatrix;
import edu.uta.futureye.algebra.FullVector;
import edu.uta.futureye.algebra.Solver;
import edu.uta.futureye.algebra.SparseBlockMatrix;
import edu.uta.futureye.algebra.SparseMatrix;
import edu.uta.futureye.algebra.SparseVector;
import edu.uta.futureye.algebra.intf.AlgebraMatrix;
import edu.uta.futureye.algebra.intf.AlgebraVector;
import edu.uta.futureye.algebra.intf.BlockMatrix;
import edu.uta.futureye.algebra.intf.Matrix;

public class TestMatrix {
	public static void testMatrixMult1() {
		SparseMatrix SA = new SparseMatrix(3,3);
		SparseMatrix SB = new SparseMatrix(3,3);
		for(int i=1;i<=3;i++) {
			for(int j=1;j<=3;j++) {
				SA.set( j,i, i+3*(j-1));
				SB.set( j,i, i+j);
			}
		}
		SA.print();
		SB.print();
		
		CompressedRowMatrix A = new CompressedRowMatrix(SA,false);
		CompressedColMatrix B = new CompressedColMatrix(SB,false);
		CompressedRowMatrix C = new CompressedRowMatrix();
		A.mult(B, C);
		C.print();
		
	}
	
	public static void testMatrixMult2() {
		SparseMatrix SA = new SparseMatrix(2,3);
		SparseMatrix SB = new SparseMatrix(3,2);
		for(int i=1;i<=3;i++) {
			for(int j=1;j<=2;j++) {
				SA.set( j,i, i+3*(j-1));
				SB.set( i,j, i+j);
			}
		}
		SA.print();
		SB.print();
		
		CompressedRowMatrix A = new CompressedRowMatrix(SA,false);
		CompressedColMatrix B = new CompressedColMatrix(SB,false);
		CompressedRowMatrix C = new CompressedRowMatrix();
		A.mult(B, C);
		C.print();
		
	}
	
	public static void testMatrixMult3() {
		SparseMatrix SA = new SparseMatrix(3,4);
		SparseMatrix SB = new SparseMatrix(4,3);
		SA.set(1, 1, 1.0);
		SA.set(2, 2, 1.0);
		SA.set(3, 3, 1.0);
		SA.set(3, 4, 1.0);
		
		SB.set(1, 1, 2.0);
		SB.set(2, 2, 2.0);
		SB.set(3, 3, 2.0);
		SB.set(4, 3, 2.0);
		
		SA.print();
		SB.print();
		
		CompressedRowMatrix A = new CompressedRowMatrix(SA,false);
		CompressedColMatrix B = new CompressedColMatrix(SB,false);
		CompressedRowMatrix C = new CompressedRowMatrix();
		A.mult(B, C);
		C.print();
		
	}
	
	public static void testMatrixMult4() {
		SparseMatrix SA = new SparseMatrix(5,8);
		SparseMatrix SB = new SparseMatrix(8,5);
		for(int i=1;i<=7;i++) {
			SA.set(1, i, i);
			SA.set(2, i, i+1.0);
		}
		SA.set(2, 8, 1.0);
		SA.set(3, 3, 2.0);
		SA.set(4, 5, 3.0);
		SA.set(4, 8, 3.0);
		SA.set(5, 1, 4.0);
		SA.set(5, 3, 4.0);
		SA.set(5, 5, 4.0);
		
		for(int i=1;i<=7;i++) {
			SB.set(i, 1, 2.0);
			SB.set(i, 5, 2.0);
		}
		SB.set(8, 5, 2.0);
		SB.set(1, 2, 2.0);
		SB.set(1, 3, 2.0);
		SB.set(3, 3, 2.0);
		SB.set(3, 4, 2.0);
		SB.set(8, 4, 2.0);
		
		SA.print();
		SB.print();
		
		CompressedRowMatrix A = new CompressedRowMatrix(SA,false);
		CompressedColMatrix B = new CompressedColMatrix(SB,false);
		CompressedRowMatrix C = new CompressedRowMatrix();
		A.mult(B, C);
		C.print();
		
		/* results:
 56.0000(1)      2.0000(2)      8.0000(3)      6.0000(4)     56.0000(5)    
 70.0000(1)      4.0000(2)     12.0000(3)     10.0000(4)     72.0000(5)    
  4.0000(1)      4.0000(3)      4.0000(4)      4.0000(5)    
  6.0000(1)      6.0000(4)     12.0000(5)    
 24.0000(1)      8.0000(2)     16.0000(3)      8.0000(4)     24.0000(5)
		 * 
		 */
		
	}
	
	
	public static void testTrans1() {
		SparseMatrix SA = new SparseMatrix(3,3);
		SparseMatrix SB = new SparseMatrix(3,3);
		for(int i=1;i<=3;i++) {
			for(int j=1;j<=3;j++) {
				SA.set( j,i, i+3*(j-1));
				SB.set( j,i, i+3*(j-1));
			}
		}
		SA.print();
		SB.print();
		
		CompressedRowMatrix A = new CompressedRowMatrix(SA,false);
		CompressedColMatrix B = new CompressedColMatrix(SB,false);
		A.getTrans().print();
		B.getTrans().print();
	}
	
	public static void testTrans2() {
		SparseMatrix SA = new SparseMatrix(5,8);
		SparseMatrix SB = new SparseMatrix(8,5);
		for(int i=1;i<=7;i++) {
			SA.set(1, i, i);
			SA.set(2, i, i+1.0);
		}
		SA.set(2, 8, 1.0);
		SA.set(3, 3, 2.0);
		SA.set(4, 5, 3.0);
		SA.set(4, 8, 3.0);
		SA.set(5, 1, 4.0);
		SA.set(5, 3, 4.0);
		SA.set(5, 5, 4.0);
		
		for(int i=1;i<=7;i++) {
			SB.set(i, 1, 2.0);
			SB.set(i, 5, 2.0);
		}
		SB.set(8, 5, 2.0);
		SB.set(1, 2, 2.0);
		SB.set(1, 3, 2.0);
		SB.set(3, 3, 2.0);
		SB.set(3, 4, 2.0);
		SB.set(8, 4, 2.0);
		
		SA.print();
		SB.print();
		
		CompressedRowMatrix A = new CompressedRowMatrix(SA,false);
		CompressedColMatrix B = new CompressedColMatrix(SB,false);
		A.getTrans().print();
		B.getTrans().print();
		
	}
	
	public static void testStorage1() {
		SparseMatrix SA = new SparseMatrix(3,3);
		SparseMatrix SB = new SparseMatrix(3,3);
		for(int i=1;i<=3;i++) {
			for(int j=1;j<=3;j++) {
				SA.set( j,i, i+3*(j-1));
				SB.set( j,i, i+3*(j-1));
			}
		}
		SA.print();
		SB.print();
		
		CompressedRowMatrix A = new CompressedRowMatrix(SA,false);
		CompressedColMatrix B = new CompressedColMatrix(SB,false);
		A.convertToCompressedCol().print();
		B.convertToCompressedRow().print();
	}
	
	public static void testStorage2() {
		SparseMatrix SA = new SparseMatrix(5,8);
		SparseMatrix SB = new SparseMatrix(8,5);
		for(int i=1;i<=7;i++) {
			SA.set(1, i, i);
			SA.set(2, i, i+1.0);
		}
		SA.set(2, 8, 1.0);
		SA.set(3, 3, 2.0);
		SA.set(4, 5, 3.0);
		SA.set(4, 8, 3.0);
		SA.set(5, 1, 4.0);
		SA.set(5, 3, 4.0);
		SA.set(5, 5, 4.0);
		
		for(int i=1;i<=7;i++) {
			SB.set(i, 1, 2.0);
			SB.set(i, 5, 2.0);
		}
		SB.set(8, 5, 2.0);
		SB.set(1, 2, 2.0);
		SB.set(1, 3, 2.0);
		SB.set(3, 3, 2.0);
		SB.set(3, 4, 2.0);
		SB.set(8, 4, 2.0);
		
		SA.print();
		SB.print();
		
		CompressedRowMatrix A = new CompressedRowMatrix(SA,false);
		CompressedColMatrix B = new CompressedColMatrix(SB,false);
		A.convertToCompressedCol().print();
		B.convertToCompressedRow().print();
		
	}

	public static void testAxpy() {
		SparseMatrix SA = new SparseMatrix(4,4);
		SparseMatrix SB = new SparseMatrix(4,4);
		SA.set(1, 1, 1.0);
		SA.set(1, 2, 2.0);
		SA.set(1, 3, 3.0);
		SA.set(1, 4, 3.0);
		SA.set(2, 1, 4.0);
		SA.set(3, 3, 4.0);
		SA.set(4, 1, 4.0);
		
		SB.set(1, 1, 2.0);
		SB.set(1, 2, 2.0);
		SB.set(1, 3, 2.0);
		SB.set(1, 4, 2.0);
		SB.set(2, 2, 2.0);
		SB.set(3, 4, 2.0);
		SB.set(4, 3, 2.0);
		
		SA.print();
		SB.print();
		
		CompressedRowMatrix A = new CompressedRowMatrix(SA,false);
		CompressedRowMatrix B = new CompressedRowMatrix(SB,false);
		A.axpy(1.0, B).print();
		
		/*
  3.0000(1c)      4.0000(2c)      5.0000(3c)      5.0000(4c)    
  4.0000(1c)      2.0000(2c)    
  4.0000(3c)      2.0000(4c)    
  4.0000(1c)      2.0000(3c)   
		 * 
		 */
		
	}	
	public static void main(String[] args) {
		//testMatrixMult1();
		//testMatrixMult2();
		//testMatrixMult3();
		//testMatrixMult4();
		//testTrans1();
		//testTrans2();
		//testStorage1();
		//testStorage2();
		testAxpy();
	}
	
	public static void testMatrixVector() {
		BlockMatrix bm = new SparseBlockMatrix(2,2);
		Matrix m11 = new SparseMatrix(3,3);
		Matrix m12 = new SparseMatrix(3,2);
		Matrix m21 = new SparseMatrix(2,3);
		Matrix m22 = new SparseMatrix(2,2);
		
		bm.setBlock(1, 1, m11);
		bm.setBlock(1, 2, m12);
		bm.setBlock(2, 1, m21);
		bm.setBlock(2, 2, m22);
		
		for(int i=1;i<=5;i++)
			for(int j=1;j<=5;j++)
				bm.set(i, j, i-1+j);
		bm.print();
		
		SparseMatrix sm = new SparseMatrix(bm.getRowDim(),bm.getColDim());
		sm.setAll(0, 0, bm.getAll());
		sm.print();
		
		SparseVector x = new SparseVector(5,3.0);
		SparseVector y = new SparseVector(5,1.0);
		SparseVector u = new SparseVector(y.getDim(),1.0);

		AlgebraMatrix aMat = new CompressedRowMatrix(sm,false);
		AlgebraVector ax = new FullVector(x);
		AlgebraVector ay = new FullVector(y);
		AlgebraVector au = new FullVector(u);
		
		Solver solver = new Solver();
		System.out.println("dot: "+x.dot(y));		
		System.out.println("dot: "+ax.dot(ay));
		
//		sm.mult(x, y);
//		y.print();
//		solver.CG1(sm, y, u);
//		u.print();
		
		aMat.mult(ax, ay);
		ay.print();
		solver.solveCG(aMat, ay, au);
		au.print();
	}
}
