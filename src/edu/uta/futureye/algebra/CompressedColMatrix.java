package edu.uta.futureye.algebra;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Map;
import java.util.Map.Entry;

import edu.uta.futureye.algebra.intf.AlgebraMatrix;
import edu.uta.futureye.algebra.intf.AlgebraVector;
import edu.uta.futureye.algebra.intf.MatrixEntry;
import edu.uta.futureye.algebra.intf.SparseMatrix;

/**
 * Compressed column matrix
 * 
 * @author liuyueming
 */
public class CompressedColMatrix implements AlgebraMatrix {
	
	/**
	 * Row indices of non-zero values. These are kept sorted within each column.
	 */
	protected int[][] rowIndex = null;
	
	protected double[][] data = null;
	
	protected int rowDim;
	protected int colDim;

	public CompressedColMatrix() {
		
	}
	
	public CompressedColMatrix(int nRow, int nCol) {
		this.colDim = nCol;
		this.rowDim = nRow;
		this.rowIndex = new int[nCol][];
		this.data = new double[nCol][];
		for(int c=0;c<this.colDim;c++) {
			this.rowIndex[c] = new int[0];
			this.data[c] = new double[0];
		}
	}
	
	public void setCol(int col,int[] rowIndex, double []data) {
		int c = col - 1;
		this.rowIndex[c] = new int[rowIndex.length];
		this.data[c] = new double[rowIndex.length];
		System.arraycopy(rowIndex, 0, this.rowIndex[c], 0, rowIndex.length);
		System.arraycopy(data, 0, this.data[c], 0, rowIndex.length);
	}
	
	public CompressedColMatrix(SparseMatrix sMat, boolean clearSparseMatrix) {
		this.rowDim = sMat.getRowDim();
		this.colDim = sMat.getColDim();
		
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
		ArrayList<IVPair>[] columns = new ArrayList[this.colDim];
		for(int i=this.colDim; --i>=0;) {
			columns[i] = new ArrayList<IVPair>();
		}
		
//		Map<Integer, Map<Integer, Double>> m = sMat.getAll();
//		for(int r=0; r<this.rowDim; r++) {
//			Map<Integer,Double> row = m.get(r+1);
//			if(row != null) {
//				for(Entry<Integer,Double> col : row.entrySet()) {
//					int c = col.getKey()-1;
//					
//					indexList[c].add(r);
//					dataList[c].add(col.getValue());
//				}
//				if(clearSparseMatrix) row.clear();
//			}
//		}
		
		for(MatrixEntry e : sMat) {
			int c = e.getCol();
			int r = e.getRow();
			double v = e.getValue();
			columns[c-1].add(new IVPair(r-1,v));
		}
		
		if(clearSparseMatrix) sMat.clearAll();
		
		this.rowIndex = new int[this.colDim][];
		this.data = new double[this.colDim][];
		for(int c=0;c<this.colDim;c++) {
			int rDim = columns[c].size();
			this.rowIndex[c] = new int[rDim];
			this.data[c] = new double[rDim];
			//sort row indices in each column
			Collections.sort(columns[c], new Comparator<IVPair>() {
				@Override
				public int compare(IVPair o1, IVPair o2) {
					if(o1.index>o2.index) return 1;
					else return -1;
				}
			});
			
			for(int r=0;r<rDim;r++) {
				IVPair pair = columns[c].get(r);
				this.rowIndex[c][r] = pair.index;
				this.data[c][r] = pair.value;
			}
			columns[c].clear();
		}
		
	}
	
	public CompressedColMatrix(SparseBlockMatrix sMat, boolean clearSparseMatrix) {
		this.rowDim = sMat.getRowDim();
		this.colDim = sMat.getColDim();
		
		Map<Integer, Map<Integer, Double>> m = sMat.getAll();
		//by column
		ArrayList<Integer>[] indexList = new ArrayList[sMat.getColDim()];
		ArrayList<Double>[] dataList = new ArrayList[sMat.getColDim()];
		for(int r=0; r<this.rowDim; r++) {
			Map<Integer,Double> row = m.get(r+1);
			if(row != null) {
				for(Entry<Integer,Double> col : row.entrySet()) {
					int c = col.getKey()-1;
					indexList[c].add(r);
					dataList[c].add(col.getValue());
				}
				if(clearSparseMatrix) row.clear();
			}
		}
		if(clearSparseMatrix) m.clear();
		
		this.rowIndex = new int[sMat.getColDim()][];
		this.data = new double[sMat.getColDim()][];
		for(int c=0;c<sMat.getColDim();c++) {
			int rDim = indexList[c].size();
			rowIndex[c] = new int[rDim];
			data[c] = new double[rDim];
			for(int r=0;r<rDim;r++) {
				rowIndex[c][r] = indexList[c].get(r);
				data[c][r] = dataList[c].get(r);
			}
			indexList[c].clear();
			dataList[c].clear();
		}
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
		throw new UnsupportedOperationException("Use CompressedRowMatrix instead for fast mult!");
	}

	@Override
	public void print() {
		for(int col=0; col<this.colDim; col++) {
			int nRow = rowIndex[col].length;
			for(int r=0; r<nRow; r++) {
				System.out.print(String.format("%8.4f(%dr)    ", 
						this.data[col][r],this.rowIndex[col][r]+1));
			}
			System.out.println();
		}
		System.out.println();
	}

	@Override
	public void mult(AlgebraMatrix B, AlgebraMatrix C) {
		throw new UnsupportedOperationException("Use CompressedRowMatrix instead for fast mult!");
	}
	
	@Override
	public AlgebraMatrix getTrans() {
		CompressedRowMatrix T = new CompressedRowMatrix();
		T.rowDim = this.colDim;
		T.colDim = this.rowDim;
		
		T.colIndex = new int[T.rowDim][];
		T.data = new double[T.rowDim][];
		
		for(int c=0; c<this.colDim; c++) {
			T.colIndex[c] = new int[this.rowIndex[c].length];
			T.data[c] = new double[this.rowIndex[c].length];
			for(int r=0; r<this.rowIndex[c].length; r++) {
				T.colIndex[c][r] = this.rowIndex[c][r];
				T.data[c][r] = this.data[c][r];
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
	 * @return A new object of class CompressedRowMatrix
	 */
	public CompressedRowMatrix getCompressedRowMatrix() {
		CompressedRowMatrix R = new CompressedRowMatrix();
		R.colDim = this.colDim;
		R.rowDim = this.rowDim;
		
		//by column
		ArrayList<Integer>[] indexList = new ArrayList[R.rowDim];
		ArrayList<Double>[] dataList = new ArrayList[R.rowDim];
		for(int i=0;i<R.rowDim;i++) {
			indexList[i] = new ArrayList<Integer>();
			dataList[i] = new ArrayList<Double>();
		}
		
		for(int c=0; c<this.colDim; c++) {
			for(int r=0; r<this.rowIndex[c].length; r++) {
				indexList[this.rowIndex[c][r]].add(c);
				dataList[this.rowIndex[c][r]].add(this.data[c][r]);
			}
		}
		
		R.colIndex = new int[R.rowDim][];
		R.data = new double[R.rowDim][];
		for(int r=0; r<R.rowDim; r++) {
			int cDim = indexList[r].size();
			R.colIndex[r] = new int[cDim];
			R.data[r] = new double[cDim];
			for(int c=0; c<cDim; c++) {
				R.colIndex[r][c] = indexList[r].get(c);
				R.data[r][c] = dataList[r].get(c);
			}
		}
		return R;
	}
	
	public void getColVector(int col, FullVector vec) {
		for(int r=0;r<this.rowDim;r++)
			vec.data[r] = 0.0;
		int c = col - 1;
		for(int r=0;r<this.rowIndex[c].length;r++)
			vec.data[this.rowIndex[c][r]] = this.data[c][r];
	}
	
	/**
	 * TODO??SparseMatrixRowMajor SparseMatrixColMajor?
	 * <p>
	 * The values are copied. So subsequent changes in this matrix 
	 * are not reflected in the returned matrix, and vice-versa.
	 * @return
	 */
	public SparseMatrix getSparseMatrix() {
		SparseMatrix rlt = new SparseMatrixRowMajor(this.rowDim,this.colDim);
		for(int c=0; c<this.colDim; c++) {
			for(int r=0; r<this.rowIndex[c].length; r++) {
				rlt.set(this.rowIndex[c][r]+1, c+1, 
						this.data[c][r]);
			}
		}
		return rlt;
	}
	
	
	/**
	 * 返回列压缩存储方式的行索引数组，行号从0开始
	 */
	public int[][] getRowIndex() {
		return this.rowIndex;
	}
	
	/**
	 * 返回列压缩存储方式的数据数组
	 */
	public double[][] getData() {
		return this.data;
	}
}
