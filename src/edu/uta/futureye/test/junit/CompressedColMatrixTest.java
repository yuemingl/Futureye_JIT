package edu.uta.futureye.test.junit;

import org.junit.Test;

import edu.uta.futureye.algebra.CompressedColMatrix;
import edu.uta.futureye.algebra.CompressedRowMatrix;
import edu.uta.futureye.algebra.SparseMatrixColMajor;
import edu.uta.futureye.algebra.SparseMatrixRowMajor;
import edu.uta.futureye.algebra.intf.SparseMatrix;

public class CompressedColMatrixTest {

	@Test
	public void testCompressedColMatrixSparseMatrixBoolean() {
		SparseMatrix SA = new SparseMatrixRowMajor(3,3);
		SparseMatrix SB = new SparseMatrixColMajor(3,3);
		for(int i=1;i<=3;i++) {
			for(int j=1;j<=3;j++) {
				SA.set( j,i, i+3*(j-1));
				SB.set( j,i, i+3*(j-1));
			}
		}
		SA.print();
		SB.print();
		
		CompressedRowMatrix A1 = new CompressedRowMatrix(SA,false);
		CompressedRowMatrix A2 = new CompressedRowMatrix(SB,false);
		CompressedColMatrix B1 = new CompressedColMatrix(SA,false);
		CompressedColMatrix B2 = new CompressedColMatrix(SB,false);
		A1.print();
		A2.print();
		B1.print();
		B2.print();
	}

}
