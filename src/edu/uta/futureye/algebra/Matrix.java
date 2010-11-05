package edu.uta.futureye.algebra;

import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;

public class Matrix {
	protected int rowDim;
	protected int colDim;
	protected Map<Integer,Map<Integer,Double>> m;
	static double eps = 1e-10;
	
	public Matrix(int rowDim, int colDim) {
		this.rowDim = rowDim;
		this.colDim = colDim;
		m = new HashMap<Integer,Map<Integer,Double>>();
	}
	
	public void set(int row, int col,double value) {
		Map<Integer,Double> arow = m.get(row);
		if(arow == null && Math.abs(value)>=eps) {
			arow = new LinkedHashMap<Integer,Double>();
			m.put(row, arow);
		} else if(arow == null && Math.abs(value)<eps) {
			return;
		}
		arow.put(col, value);
	}
	
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
	
	public void plusValue(int row, int col,double value) {
		set(row,col,get(row,col)+value);
	}
	
	public int getRowDim() {
		return rowDim;
	}
	
	public int getColDim() {
		return colDim;
	}
	
	public Map<Integer,Map<Integer,Double>> getAll() {
		return m;
	}
	
	public void print() {
		for(int i=1;i<=rowDim;i++) {
			for(int j=1;j<=colDim;j++) {
				System.out.print(String.format("%8.6f", get(i,j))+"   ");
			}
			System.out.println("");
		}
		System.out.println("");
	}
	
}
