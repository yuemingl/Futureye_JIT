package edu.uta.futureye.algebra;

import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Map.Entry;

import edu.uta.futureye.algebra.intf.BlockMatrix;
import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.util.FutureyeException;

public class SparseBlockMatrix implements BlockMatrix {
	protected int rowBlockDim;
	protected int colBlockDim;
	protected double defaultValue = 0.0;
	protected Map<Integer,Map<Integer,Matrix>> m = 
		new LinkedHashMap<Integer,Map<Integer,Matrix>>();
	protected String name="";

	/**
	 * Construct a rowBlockDim*colBlockDim block matrix
	 * @param rowBlockDim
	 * @param colBlockDim
	 */
	public SparseBlockMatrix(int rowBlockDim, int colBlockDim) {
		this.rowBlockDim = rowBlockDim;
		this.colBlockDim = colBlockDim;

	}
	
	/**
	 * Construct a rowBlockDim*colBlockDim block matrix
	 * with a default value
	 * @param rowBlockDim
	 * @param colBlockDim
	 * @param defaultValue
	 */
	public SparseBlockMatrix(int rowBlockDim, int colBlockDim,
			double defaultValue) {
		this.rowBlockDim = rowBlockDim;
		this.colBlockDim = colBlockDim;
		this.defaultValue = defaultValue;
	}	
	
	@Override
	public int getRowBlockDim() {
		return rowBlockDim;
	}

	@Override
	public int getColBlockDim() {
		return colBlockDim;
	}

	@Override
	public Matrix getBlock(int rowBlocl, int colBlock) {
		Map<Integer,Matrix> arow = m.get(rowBlocl);
		if(arow == null) {
			return null;
		} else {
			Matrix v = arow.get(colBlock);
			if(v == null) {
				return null;
			} else {
				return v;
			}
		}
	}
	
	@Override
	public void setBlock(int row, int col, Matrix subMat) {
		Map<Integer,Matrix> arow = m.get(row);
		if(arow == null) {
			arow = new LinkedHashMap<Integer,Matrix>();
			m.put(row, arow);
		}
		arow.put(col, subMat);
	}

	@Override
	public Map<Integer, Map<Integer, Matrix>> getAllBlock() {
		return m;
	}

	@Override
	public void setColDim(int nColDim) {
		throw new UnsupportedOperationException();
	}

	@Override
	public void setRowDim(int nRowDim) {
		throw new UnsupportedOperationException();
	}
	
	@Override
	public int getColDim() {
		int rlt = 0;
		Map<Integer,Matrix> row = m.get(1);
		for(Entry<Integer,Matrix> mat : row.entrySet())
			rlt += mat.getValue().getColDim();
		return rlt;
	}

	@Override
	public int getRowDim() {
		int rlt = 0;
		for(Entry<Integer,Map<Integer,Matrix>> row : m.entrySet())
			rlt += row.getValue().get(1).getRowDim();
		return rlt;
	}

	@Override
	public void set(int row, int col, double value) {
		int rBase = 0;
		Map<Integer, Matrix> rMats = null;
		boolean find = false;
		for(int r=1;r<=this.rowBlockDim;r++) {
			rMats = m.get(r);
			int rUpper = rMats.get(1).getRowDim();
			if(rBase < row && row <= rBase+rUpper) {
				find = true;
				break;
			} else {
				rBase += rUpper;
			}
		}
		if(!find) {
			throw new FutureyeException("row="+row+"; rMax="+rBase);
		}
		int cBase = 0;
		Matrix cMat = null;
		if(rMats != null) {
			find = false;
			for(int c=1;c<=this.colBlockDim;c++) {
				cMat = rMats.get(c);
				int cUpper = cMat.getColDim();
				if(cBase < col && col <= cBase+cUpper) {
					find = true;
					break;
				} else {
					cBase += cUpper;
				}
			}
			if(!find) {
				FutureyeException ex = new FutureyeException("col="+col+"; cMax="+cBase);
				ex.printStackTrace();
				System.exit(-1);
			}
		}
		if(cMat != null)
			cMat.set(row - rBase, col - cBase, value);
	}
	
	@Override
	public double get(int row, int col) {
		int rBase = 0;
		Map<Integer, Matrix> rMats = null;
		boolean find = false;
		for(int r=1;r<=this.rowBlockDim;r++) {
			rMats = m.get(r);
			int rUpper = rMats.get(1).getRowDim();
			if(rBase < row && row <= rBase+rUpper) {
				find = true;
				break;
			} else {
				rBase += rUpper;
			}
		}
		if(!find) {
			FutureyeException ex = new FutureyeException("row="+row+"; rMax="+rBase);
			ex.printStackTrace();
			System.exit(-1);
		}
		int cBase = 0;
		Matrix cMat = null;
		if(rMats != null) {
			find = false;
			for(int c=1;c<=this.colBlockDim;c++) {
				cMat = rMats.get(c);
				int cUpper = cMat.getColDim();
				if(cBase < col && col <= cBase+cUpper) {
					find = true;
					break;
				} else {
					cBase += cUpper;
				}
			}
			if(!find) {
				FutureyeException ex = new FutureyeException("col="+col+"; cMax="+cBase);
				ex.printStackTrace();
				System.exit(-1);
			}
		} else 
			return 0.0;
		
		if(cMat != null)
			return cMat.get(row - rBase, col - cBase);
		else
			return 0.0;
	}

	@Override
	public Map<Integer, Map<Integer, Double>> getAll() {
		int rBlkBase = 0;
		Matrix mat = null;
		SparseMatrix sm = new SparseMatrix(this.getRowDim(),this.getColDim());
		for(int i=1;i<=this.rowBlockDim;i++) {
			int cBlkBase = 0;
			for(int j=1;j<=this.colBlockDim;j++) {
				mat = this.getBlock(i, j);
				sm.setAll(rBlkBase,cBlkBase,mat.getAll());
				cBlkBase += mat.getColDim();
			}
			rBlkBase += mat.getRowDim();
		}
		return sm.getAll();
	}

	@Override
	public void setAll(int nRowBase, int nColBase,
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
	public void add(int row, int col, double value) {
		set(row,col,get(row,col)+value);
	}

	@Override
	public void print() {
		int rowDim = this.getRowDim();
		int colDim = this.getColDim();
		for(int i=1;i<=rowDim;i++) {
			for(int j=1;j<=colDim;j++) {
				System.out.print(String.format("%8.4f", get(i,j))+"   ");
			}
			System.out.println();
		}
		System.out.println();
	}
	
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append("SparseBlockMatrix:"+name+"(").
			append(this.rowBlockDim).
			append(",").
			append(this.colBlockDim).
			append("):N0R=").
			append(m.size()).
			append("\n");
		
		for(int i=1;i<=this.rowBlockDim;i++) {
			for(int j=1;j<=this.colBlockDim;j++) {
				sb.append("(").append(i).append(",").append(j).append(")=");
				sb.append(m.get(i).get(j)).append("\n");
			}
		}
		return sb.toString();
	}
	
	@Override
	public void mult(Vector x, Vector y) {
		throw new UnsupportedOperationException();
	}
	
	@Override
	public Matrix trans() {
		throw new UnsupportedOperationException();
	}
	
	@Override
	public Matrix copy() {
		throw new UnsupportedOperationException();
	}
	
	@Override
	public void clear() {
		throw new UnsupportedOperationException();
	}
	
	@Override
	public String getName() {
		return name;
	}

	@Override
	public Matrix setName(String name) {
		this.name = name;
		return this;
	}

}
