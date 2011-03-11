package edu.uta.futureye.algebra;

import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;

import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.Vector;

public class SparseMatrix implements Matrix {
	protected int rowDim;
	protected int colDim;
	protected double defaultValue = 0.0;
	protected Map<Integer,Map<Integer,Double>> m = 
		new HashMap<Integer,Map<Integer,Double>>();
	
	public SparseMatrix() {
	}
	
	public SparseMatrix(int rowDim, int colDim) {
		this.rowDim = rowDim;
		this.colDim = colDim;
	}
	
	public SparseMatrix(int rowDim, int colDim,double defaultValue) {
		this.rowDim = rowDim;
		this.colDim = colDim;
		this.defaultValue = defaultValue;
	}

	@Override
	public void setColDim(int nColDim) {
		this.colDim = nColDim;
	}

	@Override
	public void setRowDim(int nRowDim) {
		this.rowDim = nRowDim;
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
	public void set(int row, int col,double value) {
		Map<Integer,Double> arow = m.get(row);
		if(arow == null) {
			if(Math.abs(value) >= Matrix.zeroEps) {
				arow = new HashMap<Integer,Double>();
				m.put(row, arow);
				arow.put(col, value);
			}
		} else {
			if(Math.abs(value) < Matrix.zeroEps)
				arow.remove(col);
			else
				arow.put(col, value);
		}
	}
	
	@Override
	public double get(int row, int col) {
		Map<Integer,Double> arow = m.get(row);
		if(arow == null) {
			return 0.0;
		} else {
			Double v = arow.get(col);
			if(v == null) {
				return 0.0;
			} else {
				return v;
			}
		}
	}

	@Override
	public void add(int row, int col, double value) {
		set(row,col,get(row,col)+value);
	}
	
	@Override
	public Map<Integer, Map<Integer, Double>> getAll() {
		return m;
	}

	@Override
	public void setAll(int nRowBase,int nColBase,
			Map<Integer, Map<Integer, Double>> map) {
		for(Entry<Integer, Map<Integer, Double>> rowEentry : map.entrySet()) {
			int nRow = rowEentry.getKey();
			Map<Integer, Double> row = rowEentry.getValue();
			for(Entry<Integer, Double> entry : row.entrySet()) {
				int nCol = entry.getKey();
				set(nRowBase+nRow,nColBase+nCol,entry.getValue());
			}
		}
	}

	@Override
	public void clear() {
		this.rowDim = 0;
		this.colDim = 0;
		for(Entry<Integer,Map<Integer,Double>> row : m.entrySet()) {
			row.getValue().clear();
		}
		this.m.clear();
	}
	
	@Override
	public void mult(Vector x, Vector y) {
		for(Entry<Integer,Map<Integer,Double>> row : m.entrySet()) {
			int nRow = row.getKey();
			y.set(nRow, 0.0);
			for(Entry<Integer,Double> col : row.getValue().entrySet()) {
				int nCol = col.getKey();
				y.add(nRow, x.get(nCol)*col.getValue());
			}
		}
	}
	
	@Override
	public void print() {
		for(int i=1;i<=rowDim;i++) {
			for(int j=1;j<=colDim;j++) {
				System.out.print(String.format("%8.4f", get(i,j))+"   ");
			}
			System.out.println();
		}
		System.out.println();
	}

	////////////////////////////////////////////////////
	
	/**
	 * 返回行压缩存储方式的列索引数组，列号从0开始
	 */
	public int [][] getColIndex() {
		int[][] colIndex = new int[m.size()][];
		for(int r=0; r<this.rowDim; r++) {
			Map<Integer,Double> row = m.get(r+1);
			colIndex[r] = new int[row.size()];
			int c = 0;
			for(Entry<Integer,Double> col : row.entrySet()) {
				int nCol = col.getKey();
				colIndex[r][c] = nCol-1;
				c++;
			}
		}
		return colIndex;
	}
}
