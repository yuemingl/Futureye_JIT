package edu.uta.futureye.algebra;

import java.util.Map;
import java.util.Map.Entry;

import edu.uta.futureye.algebra.intf.AlgebraMatrix;
import edu.uta.futureye.algebra.intf.AlgebraVector;

public class CompressedRowMatrix implements AlgebraMatrix {
	protected int[][] colIndex = null;
	protected double[][] data = null;
	protected int rowDim;
	protected int colDim;

	public CompressedRowMatrix(SparseMatrix sMat, boolean clearSparseMatrix) {
		this.rowDim = sMat.getRowDim();
		this.colDim = sMat.getColDim();
		
		Map<Integer, Map<Integer, Double>> m = sMat.getAll();
		this.colIndex = new int[m.size()][];
		this.data = new double[m.size()][];
		for(int r=0; r<this.rowDim; r++) {
			Map<Integer,Double> row = m.get(r+1);
			this.colIndex[r] = new int[row.size()];
			this.data[r] = new double[row.size()];
			int c = 0;
			for(Entry<Integer,Double> col : row.entrySet()) {
				this.colIndex[r][c] = col.getKey()-1;
				this.data[r][c] = col.getValue();
				c++;
			}
			if(clearSparseMatrix) row.clear();
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
				System.out.print(String.format("%8.4f(%i)", 
						this.data[row][c],this.colIndex[row][c])+"   ");
			}
			System.out.println();
		}
		System.out.println();
	}
}
