/**
 * Copyright (c) 2010, nkliuyueming@gmail.com. All rights reserved.
 * 
 * 
 */
package edu.uta.futureye.algebra;

import edu.uta.futureye.algebra.intf.AlgebraMatrix;
import edu.uta.futureye.algebra.intf.AlgebraVector;
import edu.uta.futureye.algebra.intf.MatrixEntry;
import edu.uta.futureye.algebra.intf.SparseMatrix;

/**
 * Full(dense) matrix, row-major storage.
 * 
 * <p>
 * TODO
 * 是否需要column-major storage的full matrix实现？
 * 
 * @author liuyueming
 *
 */
public class FullMatrix implements AlgebraMatrix {
	
	/**
	 * Number of rows
	 */
	protected int rowDim;
	
	/**
	 * Number of columns
	 */
	protected int colDim;	
	
	/**
	 * <tt>data</tt> has a row-major formatting:<br>
     *  <br>
     * data[ row ][ column ]
     * <p>
     * Row and column indices are all from 0
	 */
	protected double[][] data = null;


	/**
     * Construct a matrix with the values and dimensions defined by the 2D array <tt>data</tt>.
     * It is assumed that <tt>data</tt> has a row-major formatting:<br>
     *  <br>
     * data[ row ][ column ]
     * 
     * @param data 2D array representation of the matrix.
	 * @param bCopy Indicates whether to copy the data or not.
	 */
	public FullMatrix(double[][] data, boolean bCopy) {
		this.rowDim = data.length;
		this.colDim = data[0].length;
		if(bCopy) {
			int nRow = this.rowDim;
			int nCol = this.colDim;
			this.data = new double[nRow][nCol];
			for(int i=0;i<nRow;i++) {
				System.arraycopy(data[i], 0, this.data[i], 0, nCol);
			}
		} else {
			this.data = data;
		}
	}
	
	/**
     * Construct a <tt>nRow*nCol</tt> matrix with zero values
     * 
	 * @param nRow Number of rows
	 * @param nCol Number of columns
	 */
	public FullMatrix(int nRow, int nCol) {
		this.rowDim = nRow;
		this.colDim = nCol;
		
		this.data = new double[nRow][];
		for(int r=0; r<nRow; r++) {
			this.data[r] = new double[nCol];
			for(int i=0;i<nCol;i++) {
				data[r][i] = 0.0;
			}
		}
	}
	
	/**
     * Construct a <tt>nRow*nCol</tt> matrix with the values defined by the 1D array <tt>dataArray</tt>.
     * It is assumed that <tt>dataArray</tt> has a row-major formatting. The data is copied.
     * 
	 * @param nRow Number of rows
	 * @param nCol Number of columns
	 * @param dataArray 1D array representation of the matrix.
	 */
	public FullMatrix(int nRow, int nCol, double[] dataArray) {
		this.rowDim = nRow;
		this.colDim = nCol;
		
		this.data = new double[nRow][];
		for(int r=0; r<nRow; r++) {
			this.data[r] = new double[nCol];
			System.arraycopy(dataArray, r*nCol, this.data[r], 0, nCol);
		}
	}
	

	/**
	 * Construct a matrix with a SparseMatrix. 
	 * The values are copied. So subsequent changes in the SparseMatrix 
	 * are not reflected in the matrix, and vice-versa.
	 * 
	 * @param sMat
	 */
	public FullMatrix(SparseMatrix sMat) {
		int m = sMat.getRowDim();
		int n = sMat.getColDim();
		this.rowDim = m;
		this.colDim = n;
		
		this.data = new double[m][n];
		for(int i=m;--i>=0;) 
			for(int j=n; --j>=0;)
				this.data[i][j] = 0.0;
		for(MatrixEntry e : sMat) {
			this.data[e.getRow()-1][e.getCol()-1] = e.getValue();
		}
//		Map<Integer, Map<Integer, Double>> m = sMat.getAll();
//		for(int r=0; r<this.rowDim; r++) {
//			this.data[r] = new double[this.colDim];
//			for(int i=0;i<this.colDim;i++) {
//				data[r][i] = 0.0;
//			}
//			Map<Integer, Double> row = m.get(r+1);
//			if(row != null) {
//				for(Entry<Integer,Double> e : row.entrySet()) {
//					this.data[r][e.getKey()-1] = e.getValue();
//				}
//			}
//		}
	}
	
	@Override
	public int getRowDim() {
		return rowDim;
	}

	@Override
	public int getColDim() {
		return colDim;
	}

	@Override
	public void mult(AlgebraVector x, AlgebraVector y) {
		double[] xData = x.getData();
		double[] yData = y.getData();
		double sum;
		for(int i=this.rowDim; --i>=0;) {
			sum = 0.0;
			for(int j=this.colDim; --j>=0;)
				sum += this.data[i][j]*xData[j];//9/29/2012 bugfix: = => +=
			yData[i] = sum;
		}
	}

	@Override
	public void mult(AlgebraMatrix B, AlgebraMatrix C) {
		if(B instanceof FullMatrix && C instanceof FullMatrix) {
			FullMatrix BB = (FullMatrix)B;
			FullMatrix CC = (FullMatrix)C;
			for(int i=0;i<this.rowDim;i++) {
				double[] pRowAA = this.data[i];
				double[] pRowCC = CC.data[i];
				for(int j=0;j<BB.colDim;j++) {
					pRowCC[j] = 0.0;
					for(int k=0;k<this.rowDim;k++) {
						pRowCC[j] += pRowAA[k]*BB.data[k][j];
					}
				}
			}
		}
	}

	/**
	 * A=A'
	 * <p>
	 * The values are copied. So subsequent changes in returned matrix are not reflected in the matrix, and vice-versa.
	 * 
	 * @return A new object of <tt>FullMatrix<tt> containing A'
	 */
	@Override
	public AlgebraMatrix getTrans() {
		FullMatrix rlt = new FullMatrix(data,true);
		for(int i=this.rowDim; --i>=0;)
			for(int j=this.colDim; --j>=0;)
				rlt.data[j][i] = this.data[i][j];
		return rlt;
	}

	@Override
	public void print() {
		int nRow = this.rowDim;
		int nCol = this.colDim;
		for(int i=0;i<nRow;i++) {
			for(int j=0;j<nCol;j++) {
				System.out.print(String.format("%8.6f   ", this.data[i][j]));
			}
			System.out.println();
		}
		System.out.println();

	}
	
	/**
	 * Deep copy of this matrix
	 * <p>
	 * The values are copied. So subsequent changes in returned matrix are not reflected in the matrix, and vice-versa.
	 * 
	 * @return A new object of <tt>FullMatrix<tt> containing the same values and dimensions as this matrix
	 */
	public FullMatrix copy() {
		FullMatrix rlt = new FullMatrix(data,true);
		return rlt;
	}

	/**
	 * Return a 2-dimensional array containing the values of matrix, the returned values of array has the form:
	 * <p>
	 * <b>values[row][column]</b>
	 * <p>
	 * The returned values is backed by this matrix, so changes in the returned array are reflected in this matrix, and vice-versa.
	 * @return Values of matrix in a 2-dimensional array
	 */
	public double[][] getData() {
		return data;
	}

}
