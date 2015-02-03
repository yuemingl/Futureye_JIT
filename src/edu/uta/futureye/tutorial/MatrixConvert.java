package edu.uta.futureye.tutorial;

import org.ejml.data.DenseMatrix64F;

import no.uib.cipr.matrix.sparse.CompColMatrix;
import no.uib.cipr.matrix.sparse.CompRowMatrix;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.impl.SparseDoubleMatrix2D;
import edu.uta.futureye.algebra.CompressedColMatrix;
import edu.uta.futureye.algebra.CompressedRowMatrix;
import edu.uta.futureye.algebra.FullMatrix;
import edu.uta.futureye.algebra.SparseMatrixRowMajor;
import edu.uta.futureye.algebra.intf.MatrixEntry;
import edu.uta.futureye.algebra.intf.SparseMatrix;

public class MatrixConvert {
	
	//----------------------Jama------------------------------
	public static Jama.Matrix toJamaMatrix(SparseMatrix A) {
		FullMatrix tmpA = new FullMatrix(A);
		Jama.Matrix A2 = new Jama.Matrix(
				tmpA.getData(),A.getRowDim(),A.getColDim());
		return A2;
	}
	
	//----------------------Colt--------------------------------
	public static DenseDoubleMatrix2D toColtDenseDoubleMatrix2D(SparseMatrix A) {
		FullMatrix tmpA = new FullMatrix(A);
		DenseDoubleMatrix2D A2 = new DenseDoubleMatrix2D(tmpA.getData());
		return A2;
	}	
	
	public static SparseDoubleMatrix2D toColtSparseDoubleMatrix2D(SparseMatrix A) {
		SparseDoubleMatrix2D A2 = new SparseDoubleMatrix2D(
				A.getRowDim(),A.getColDim());
		for(MatrixEntry e : A) {
			A2.setQuick(e.getRow(), e.getCol(), e.getValue());
		}
		return A2;
	}	
	
	
	//---------------------------MTJ----------------------------
	
	public static CompRowMatrix toMTJCompRowMatrix(SparseMatrix A) {
        CompressedRowMatrix tmpA = new CompressedRowMatrix(A,false);
        CompRowMatrix cnvtA = new CompRowMatrix(
        		A.getRowDim(),A.getColDim(),tmpA.getColIndex());
        for(MatrixEntry e : A) {
            cnvtA.set(e.getRow()-1, e.getCol()-1, e.getValue());
        }
        return cnvtA;
	}
	
	public static CompColMatrix toMTJCompColMatrix(SparseMatrix A) {
        CompressedColMatrix tmpA = new CompressedColMatrix(A,false);
        CompColMatrix cnvtA = new CompColMatrix(
        		A.getRowDim(),A.getColDim(),tmpA.getRowIndex());
        for(MatrixEntry e : A) {
            cnvtA.set(e.getRow()-1, e.getCol()-1, e.getValue());
        }
        return cnvtA;
	}
	
	public static SparseMatrix fromMTJ(CompRowMatrix A) {
		SparseMatrix A2 = new SparseMatrixRowMajor(
				A.numRows(),A.numColumns());
        for(no.uib.cipr.matrix.MatrixEntry e : A) {
            A2.set(e.row()+1, e.column()+1, e.get());
        }
        return A2;
	}
	
	public static SparseMatrix fromMTJ(CompColMatrix A) {
		SparseMatrix A2 = new SparseMatrixRowMajor(
				A.numRows(),A.numColumns());
        for(no.uib.cipr.matrix.MatrixEntry e : A) {
            A2.set(e.row()+1, e.column()+1, e.get());
        }
        return A2;
	}

	
	//---------------------------EJML------------------------
	public static DenseMatrix64F toEJMLDenseMatrix64F(SparseMatrix A) {
		FullMatrix tmpA = new FullMatrix(A);
		DenseMatrix64F A2 = new DenseMatrix64F(tmpA.getData());
		return A2;
	}	
	

	/**
	 * @param args
	 */
	public static void main(String[] args) {
	}
	
}
