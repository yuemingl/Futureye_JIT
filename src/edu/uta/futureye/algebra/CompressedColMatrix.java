package edu.uta.futureye.algebra;

import java.util.ArrayList;
import java.util.Map;
import java.util.Map.Entry;

import edu.uta.futureye.algebra.intf.AlgebraMatrix;
import edu.uta.futureye.algebra.intf.AlgebraVector;

public class CompressedColMatrix implements AlgebraMatrix {
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
	}
	
	public void setCol(int col,int[] rowIndex, double []data) {
		int c = col - 1;
		this.rowIndex[c] = new int[rowIndex.length];
		this.data[c] = new double[rowIndex.length];
		for(int r=0; r<rowIndex.length; r++) {
			this.rowIndex[c][r] = rowIndex[r];
			this.data[c][r] = data[r];
		}
	}
	
	public CompressedColMatrix(SparseMatrix sMat, boolean clearSparseMatrix) {
		this.rowDim = sMat.getRowDim();
		this.colDim = sMat.getColDim();
		
		Map<Integer, Map<Integer, Double>> m = sMat.getAll();
		//by column
		ArrayList<Integer>[] indexList = new ArrayList[sMat.colDim];
		ArrayList<Double>[] dataList = new ArrayList[sMat.colDim];
		for(int i=0;i<sMat.colDim;i++) {
			indexList[i] = new ArrayList<Integer>();
			dataList[i] = new ArrayList<Double>();
		}
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
		
		this.rowIndex = new int[sMat.colDim][];
		this.data = new double[sMat.colDim][];
		for(int c=0;c<sMat.colDim;c++) {
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
		throw new UnsupportedOperationException();
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
		throw new UnsupportedOperationException();

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
	
	public CompressedRowMatrix convertToCompressedRow() {
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
	
	public SparseMatrix getSparseMatrix() {
		SparseMatrix rlt = new SparseMatrix();
		for(int c=0; c<this.colDim; c++) {
			for(int r=0; r<this.rowIndex[c].length; r++) {
				rlt.set(this.rowIndex[c][r], c+1, 
						this.data[c][r]);
			}
		}
		return rlt;
	}
}
