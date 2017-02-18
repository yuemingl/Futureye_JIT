/**
 * Copyright (c) 2010, nkliuyueming@gmail.com. All rights reserved.
 * 
 * 
 */
package edu.uta.futureye.algebra;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;

import edu.uta.futureye.algebra.intf.AlgebraMatrix;
import edu.uta.futureye.algebra.intf.AlgebraVector;
import edu.uta.futureye.algebra.intf.MatrixEntry;
import edu.uta.futureye.algebra.intf.SparseMatrix;
import edu.uta.futureye.util.Constant;

/**
 * Compressed row matrix
 * 
 * @author liuyueming
 */
public class CompressedRowMatrix implements AlgebraMatrix {
	
	/**
	 * Column indices of non-zero values. These are kept sorted within each row.
	 */
	protected int[][] colIndex = null;

	protected double[][] data = null;
	
	protected int rowDim;
	protected int colDim;

	public CompressedRowMatrix() {
		
	}
	
	public CompressedRowMatrix(int nRow, int nCol) {
		this.colDim = nCol;
		this.rowDim = nRow;
		this.colIndex = new int[nRow][];
		this.data = new double[nRow][];
		for(int r=0;r<this.rowDim;r++) {
			this.colIndex[r] = new int[0];
			this.data[r] = new double[0];
		}
	}
	
	public void setRow(int row,int[] colIndex, double []data) {
		int r = row - 1;
		this.colIndex[r] = new int[colIndex.length];
		this.data[r] = new double[colIndex.length];
		System.arraycopy(colIndex, 0, this.colIndex[r], 0, colIndex.length);
		System.arraycopy(data, 0, this.data[r], 0, colIndex.length);
	}
	
	public CompressedRowMatrix(SparseMatrix sMat, boolean clearSparseMatrix) {
		this.rowDim = sMat.getRowDim();
		this.colDim = sMat.getColDim();
		
//		Map<Integer, Map<Integer, Double>> m = sMat.getAll();
//		this.colIndex = new int[this.rowDim][];
//		this.data = new double[this.rowDim][];
//		for(int r=0; r<this.rowDim; r++) {
//			Map<Integer,Double> row = m.get(r+1);
//			int size = 0;
//			if(row != null) size = row.size();
//			this.colIndex[r] = new int[size];
//			this.data[r] = new double[size];
//			int c = 0;
//			if(row != null) {
//				for(Entry<Integer,Double> col : row.entrySet()) {
//					this.colIndex[r][c] = col.getKey()-1;
//					this.data[r][c] = col.getValue();
//					c++;
//				}
//			}
//			if(clearSparseMatrix) if(row!=null) row.clear();
//		}
//		if(clearSparseMatrix) m.clear();
		
		class IVPair {
			int index;
			double value;
			IVPair(int index, double value) {
				this.index = index;
				this.value = value;
			}
		}
		
		//by column
		@SuppressWarnings("unchecked")
		ArrayList<IVPair>[] rows = new ArrayList[this.rowDim];
		for(int i=this.rowDim; --i>=0;) {
			rows[i] = new ArrayList<IVPair>();
		}
		for(MatrixEntry e : sMat) {
			int c = e.getCol();
			int r = e.getRow();
			double v = e.getValue();
			rows[r-1].add(new IVPair(c-1,v));
		}
		
		if(clearSparseMatrix) sMat.clearAll();
		
		this.colIndex = new int[this.colDim][];
		this.data = new double[this.colDim][];
		for(int r=this.rowDim; --r>=0;) {
			int cDim = rows[r].size();
			this.colIndex[r] = new int[cDim];
			this.data[r] = new double[cDim];
			//sort column indices in each row
			Collections.sort(rows[r], new Comparator<IVPair>() {
				@Override
				public int compare(IVPair o1, IVPair o2) {
					if(o1.index>o2.index) return 1;
					else return -1;
				}
			});
			
			for(int c=cDim; --c>=0;) {
				IVPair pair = rows[r].get(c);
				this.colIndex[r][c] = pair.index;
				this.data[r][c] = pair.value;
			}
			rows[r].clear();
		}		
	}
	
	public CompressedRowMatrix(SparseBlockMatrix sMat, boolean clearSparseMatrix) {
		this.rowDim = sMat.getRowDim();
		this.colDim = sMat.getColDim();
		
		Map<Integer, Map<Integer, Double>> m = sMat.getAll();
		this.colIndex = new int[this.rowDim][];
		this.data = new double[this.rowDim][];
		for(int r=0; r<this.rowDim; r++) {
			Map<Integer,Double> row = m.get(r+1);
			int size = 0;
			if(row != null) size = row.size();
			this.colIndex[r] = new int[size];
			this.data[r] = new double[size];
			int c = 0;
			if(row != null) {
				for(Entry<Integer,Double> col : row.entrySet()) {
					this.colIndex[r][c] = col.getKey()-1;
					this.data[r][c] = col.getValue();
					c++;
				}
			}
			if(clearSparseMatrix) if(row!=null) row.clear();
		}
		if(clearSparseMatrix) m.clear();
	}
	
	@Override
	public int getColDim() {
		return this.colDim;
	}

	@Override
	public int getRowDim() {
		return this.rowDim;
	}

	@Override
	public void mult(AlgebraVector x, AlgebraVector y) {
		double[] xData = x.getData();
		double[] yData = y.getData();
		for(int row=0; row<this.rowDim; row++) {
			int nCol = colIndex[row].length;
			double v = 0.0;
			for(int c=0; c<nCol; c++) {
				v += data[row][c] * xData[colIndex[row][c]];
			}
			yData[row] = v;
		}
	}

	@Override
	public void print() {
		for(int row=0; row<this.rowDim; row++) {
			int nCol = colIndex[row].length;
			for(int c=0; c<nCol; c++) {
				System.out.print(String.format("%8.4f(%dc)    ", 
						this.data[row][c],this.colIndex[row][c]+1));
			}
			System.out.println();
		}
		System.out.println();
	}
	

	@Override
	public void mult(AlgebraMatrix B, AlgebraMatrix C) {
		if(B instanceof CompressedColMatrix && C instanceof CompressedRowMatrix) {
			double[] tmpData = new double[this.colDim]; //保存一整行数据的临时数组
			int[] tmpFlag = new int[this.colDim]; //保存一整行标记的临时数组
			for(int j=0;j<this.colDim;j++) tmpData[j] = 0.0;
			for(int j=0;j<this.colDim;j++) tmpFlag[j] = -1;
			
			CompressedColMatrix BB = (CompressedColMatrix)B;
			CompressedRowMatrix CC = (CompressedRowMatrix)C;
			CC.rowDim = this.rowDim;
			CC.colDim = BB.colDim;
			CC.colIndex = new int[CC.rowDim][];
			CC.data = new double[CC.rowDim][];
	
			int[] cColIndex= new int[BB.colDim];
			double[] cColData = new double[BB.colDim];
			for(int row=0; row<this.rowDim; row++) {
				
				int nCol = this.colIndex[row].length;
				for(int c=0; c<nCol; c++) {
					int ci = this.colIndex[row][c];
					tmpData[ci] = this.data[row][c]; //非零列的数据保存在临时数组中
					tmpFlag[ci] = row; //非零列的列号标记
				}
				int total = 0; //非零列总数计数器置0
				for(int col=0; col<BB.colDim; col++) {
					double v = 0.0;	
					int nRow = BB.rowIndex[col].length;
					double[] pCol = BB.data[col];
					int[] pColIdx = BB.rowIndex[col];
					for(int r=0; r<nRow; r++) {
						int ri = pColIdx[r];
						if(tmpFlag[ri] == row)
							v += tmpData[ri]*pCol[r];
					}
					if(Math.abs(v) > Constant.eps) {
						cColData[total] = v;
						cColIndex[total] = col;
						total++;
					}
				}
				CC.colIndex[row] = new int[total];
				CC.data[row] = new double[total];
				System.arraycopy(cColIndex, 0, CC.colIndex[row], 0, total);
				System.arraycopy(cColData, 0, CC.data[row], 0, total);
			}
		} else if(B instanceof FullMatrix && C instanceof CompressedRowMatrix) {
			FullMatrix BB = (FullMatrix)B;
			CompressedRowMatrix CC = (CompressedRowMatrix)C;
			CC.rowDim = this.rowDim;
			CC.colDim = BB.colDim;
			CC.colIndex = new int[CC.rowDim][];
			CC.data = new double[CC.rowDim][];
			
			double[] cColData = new double[BB.colDim];
			int[] cColIndex = new int[BB.colDim];
			
			for(int row=0; row<this.rowDim; row++) {
				int nCol = this.colIndex[row].length;
				double[] pRow = this.data[row];
				int[] pRowIdx = this.colIndex[row];
				int total = 0;
				for(int col=0; col<BB.colDim; col++) {
					double v = 0.0;	
					for(int c=0; c<nCol; c++) {
						int ci = pRowIdx[c];
						v += pRow[c]*BB.data[ci][col];
					}
					if(Math.abs(v) > Constant.eps) {
						cColData[total] = v;
						cColIndex[total] = col;
						total++;
					}
				}
				CC.colIndex[row] = new int[total];
				CC.data[row] = new double[total];
				System.arraycopy(cColIndex, 0, CC.colIndex[row], 0, total);
				System.arraycopy(cColData, 0, CC.data[row], 0, total);
			}
		} else {
			throw new UnsupportedOperationException();
		}
	}

	@Override
	public AlgebraMatrix getTrans() {
		CompressedColMatrix T = new CompressedColMatrix();
		T.rowDim = this.colDim;
		T.colDim = this.rowDim;
		
		T.rowIndex = new int[T.colDim][];
		T.data = new double[T.colDim][];
		
		for(int r=0; r<this.rowDim; r++) {
			T.rowIndex[r] = new int[this.colIndex[r].length];
			T.data[r] = new double[this.colIndex[r].length];
			for(int c=0; c<this.colIndex[r].length; c++) {
				T.rowIndex[r][c] = this.colIndex[r][c];
				T.data[r][c] = this.data[r][c];
			}
		}
		return T;		
	}

	/**
	 * CompressedRowMatrix与CompressedColMatrix视为一对压缩存储的矩阵类，
	 * 他们之间可以通过该函数互相转换
	 * <p>
	 * The values are copied. So subsequent changes in this matrix 
	 * are not reflected in the returned matrix, and vice-versa.
	 * 
	 * @return A new object of class CompressedRowMatrix, data are copied.
	 */
	public CompressedColMatrix getCompressedColMatrix() {
		CompressedColMatrix C = new CompressedColMatrix();
		C.colDim = this.colDim;
		C.rowDim = this.rowDim;
		
		//by column
		ArrayList<Integer>[] indexList = new ArrayList[C.colDim];
		ArrayList<Double>[] dataList = new ArrayList[C.colDim];
		for(int i=0;i<C.colDim;i++) {
			indexList[i] = new ArrayList<Integer>();
			dataList[i] = new ArrayList<Double>();
		}
		
		for(int r=0; r<this.rowDim; r++) {
			for(int c=0; c<this.colIndex[r].length; c++) {
				indexList[this.colIndex[r][c]].add(r);
				dataList[this.colIndex[r][c]].add(this.data[r][c]);
			}
		}
		
		C.rowIndex = new int[C.colDim][];
		C.data = new double[C.colDim][];
		for(int c=0; c<C.colDim; c++) {
			int rDim = indexList[c].size();
			C.rowIndex[c] = new int[rDim];
			C.data[c] = new double[rDim];
			for(int r=0; r<rDim; r++) {
				C.rowIndex[c][r] = indexList[c].get(r);
				C.data[c][r] = dataList[c].get(r);
			}
		}
		return C;
	}
	
	public void getRowVector(int row, FullVector vec) {
		for(int r=0;r<this.colDim;r++)
			vec.data[r] = 0.0;
		int r = row - 1;
		int[] idx = this.colIndex[r];
		double[] dat = this.data[r];
		for(int c=0; c<idx.length; c++)
			vec.data[idx[c]] = dat[c];
	}
	
	/**
	 * A = a*A+B
	 * @param a
	 * @param B
	 * @return
	 */
	public CompressedRowMatrix axpy(double a, CompressedRowMatrix B) {
		double[] tmpData = new double[this.colDim];
		int[] tmpFlag = new int[this.colDim];
		for(int j=0;j<this.colDim;j++) tmpData[j] = 0.0;
		for(int j=0;j<this.colDim;j++) tmpFlag[j] = -1;
		
		Map<Integer, Double> mapRow = new HashMap<Integer, Double>();
		for(int row=0; row<this.rowDim; row++) {
			mapRow.clear();
			
			int nCol = this.colIndex[row].length;
			for(int c=0; c<nCol; c++) {
				int ci = this.colIndex[row][c];
				Double v = mapRow.get(ci);
				if(v == null) {
					mapRow.put(ci, a*this.data[row][c]);
				} else {
					mapRow.put(ci, v + a*this.data[row][c]);
				}
			}
			
			nCol = B.colIndex[row].length;
			for(int c=0; c<nCol; c++) {
				int ci = B.colIndex[row][c];
				Double v = mapRow.get(ci);
				if(v == null) {
					mapRow.put(ci, B.data[row][c]);
				} else {
					mapRow.put(ci, v + B.data[row][c]);
				}
			}
			this.colIndex[row] = new int[mapRow.size()];
			this.data[row] = new double[mapRow.size()];
			int c = 0;
			for(Entry<Integer,Double> ety : mapRow.entrySet()) {
				this.colIndex[row][c] = ety.getKey();
				this.data[row][c] = ety.getValue();
				c++;
			}
		}
		return this;
	}
	
	
	public CompressedRowMatrix ax(double a) {
		for(int row=0; row<this.rowDim; row++) {
			for(int c=0; c<this.colIndex[row].length; c++)
				this.data[row][c] *= a;
		}
		return this;
	}
	
	/**
	 * 
	 * <p>
	 * The values are copied. So subsequent changes in this matrix 
	 * are not reflected in the returned matrix, and vice-versa.
	 * 
	 * @return
	 */
	public SparseMatrix getSparseMatrix() {
		SparseMatrix rlt = new SparseMatrixRowMajor(this.rowDim,this.colDim);
		for(int r=0; r<this.rowDim; r++) {
			for(int c=0; c<this.colIndex[r].length; c++) {
				rlt.set(r+1, this.colIndex[r][c]+1,
						this.data[r][c]);
			}
		}
		return rlt;
	}
	
	/**
	 * 返回行压缩存储方式的列索引数组，列号从0开始
	 */
	public int[][] getColIndex() {
		return this.colIndex;
	}
	
	/**
	 * 返回行压缩存储方式的数据数组
	 */
	public double[][] getData() {
		return this.data;
	}

}
