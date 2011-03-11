package edu.uta.futureye.test;

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
import edu.uta.futureye.algebra.intf.Vector;

public class TestMatrix {
	public static void main(String[] args) {
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
		solver.CG(aMat, ay, au);
		au.print();
	}
}
