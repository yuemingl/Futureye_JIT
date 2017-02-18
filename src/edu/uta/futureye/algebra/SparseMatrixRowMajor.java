/**
 * Copyright (c) 2010, nkliuyueming@gmail.com. All rights reserved.
 * 
 * 
 */
package edu.uta.futureye.algebra;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Map.Entry;

import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.MatrixEntry;
import edu.uta.futureye.algebra.intf.SparseMatrix;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.io.MatlabMatFileWriter;
import edu.uta.futureye.util.FutureyeException;
import edu.uta.futureye.util.Sequence;

/**
 * Row-major map based storage sparse matrix implementation
 * <p>
 * This implementation provides constant-time performance for the basic operations (get and set), 
 * assuming the background hash function disperses the elements properly among the buckets.
 * 
 * @author liuyueming
 *
 */
public class SparseMatrixRowMajor implements SparseMatrix {
	protected int rowDim = 0;
	protected int colDim = 0;
	protected double defaultValue = 0.0;
	
	protected Map<Integer,Map<Integer,Double>> m = 
		new HashMap<Integer,Map<Integer,Double>>();
	
	protected String name = this.getClass().getSimpleName()+Sequence.getInstance().nextSeq();
	
	public SparseMatrixRowMajor() {
	}
	
	public SparseMatrixRowMajor(String name) {
		this.name = name;
	}
	
	public SparseMatrixRowMajor(int rowDim, int colDim) {
		this.rowDim = rowDim;
		this.colDim = colDim;
	}
	
	public SparseMatrixRowMajor(String name, int rowDim, int colDim) {
		this.name = name;
		this.rowDim = rowDim;
		this.colDim = colDim;
	}
	
	public SparseMatrixRowMajor(int rowDim, int colDim,double defaultValue) {
		this.rowDim = rowDim;
		this.colDim = colDim;
		this.defaultValue = defaultValue;
	}
	
	public SparseMatrixRowMajor(String name,int rowDim, int colDim,double defaultValue) {
		this.name = name;
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
		if(rowDim != 0) {
			if(row < 1 || row > rowDim)
				throw new FutureyeException("Row number "+row+" exceeds dimenstion [1,"+rowDim+"]");
		}
		if(colDim != 0) {
			if(col < 1 || col > colDim)
				throw new FutureyeException("Column number "+col+" exceeds dimenstion [1,"+colDim+"]");
		}
		Map<Integer,Double> aRow = m.get(row);
		if(aRow == null) {
			if(Math.abs(value) >= Matrix.zeroEps) {
				aRow = new HashMap<Integer,Double>();
				m.put(row, aRow);
				aRow.put(col, value);
			}
		} else {
			if(Math.abs(value) < Matrix.zeroEps)
				aRow.remove(col);
			else
				aRow.put(col, value);
		}
	}
	
	@Override
	public double get(int row, int col) {
		if(rowDim != 0) {
			if(row < 1 || row > rowDim)
				throw new FutureyeException("Row number "+row+" exceeds dimenstion [1,"+rowDim+"]");
		}
		if(colDim != 0) {
			if(col < 1 || col > colDim)
				throw new FutureyeException("Column number "+col+" exceeds dimenstion [1,"+colDim+"]");
		}
		Map<Integer,Double> aRow = m.get(row);
		if(aRow == null) {
			return 0.0;
		} else {
			Double v = aRow.get(col);
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
	public void clearAll() {
		this.rowDim = 0;
		this.colDim = 0;
		this.defaultValue = 0.0;
		//this.name = null;
		for(Entry<Integer,Map<Integer,Double>> row : m.entrySet()) {
			row.getValue().clear();
		}
		this.m.clear();
	}
	
	@Override
	public void clearData() {
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
	
	/**
	 * A=A' 
	 * 
	 */
	@Override
	public SparseMatrixRowMajor trans() {
		Map<Integer,Map<Integer,Double>> m2 = m;
		m = new HashMap<Integer,Map<Integer,Double>>();
		int dim = this.colDim;
		this.colDim = this.rowDim;
		this.rowDim = dim;
		for(Entry<Integer,Map<Integer,Double>> row : m2.entrySet()) {
			int nRow = row.getKey();
			for(Entry<Integer,Double> col : row.getValue().entrySet()) {
				int nCol = col.getKey();
				set(nCol, nRow, col.getValue());
			}
		}
		return this;
	}
	
	/**
	 * An overriding method can also return a subtype of the type returned by the overridden method. 
	 * This is called a covariant return type.
	 */
	@Override
	public SparseMatrixRowMajor copy() {
		SparseMatrixRowMajor newM = new SparseMatrixRowMajor(this.rowDim,this.colDim);
		newM.setAll(0, 0, this.m);
		return newM;
	}
	
	@Override
	public void print() {
		for(int i=1;i<=rowDim;i++) {
			for(int j=1;j<=colDim;j++) {
				System.out.print(String.format("%8.6f   ", get(i,j)));
			}
			System.out.println();
		}
		System.out.println();
	}
	
	public String toString() {
		return "SparseMatrix:"+name+"("+
			this.rowDim+","+this.colDim+
			"):N0R="+m.size();
	}

	////////////////////////////////////////////////////
	
	/**
	 * 返回行压缩存储方式的列索引数组，列号从0开始
	 */
	public int [][] getColIndex() {
		int[][] colIndex = new int[m.size()][];
		//for(int r=0; r<this.rowDim; r++) {
		for(int r=this.rowDim; --r>=0;) {
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
	

	@Override
	public String getName() {
		return name;
	}

	@Override
	public SparseMatrix setName(String name) {
		this.name = name;
		return this; 
	}

	/**
	 * Get number of non zero values
	 * 
	 */
	public int getNonZeroNumber() {
		int rlt = 0;
		for(Entry<Integer,Map<Integer,Double>> e1 : m.entrySet()) {
			rlt += e1.getValue().size();
		}
		return rlt;
	}
	
	/**
	 * Return the product A[row,:]*B[:,col]
	 * where <tt>A==this</tt>
	 * 
	 * @param B
	 * @param row
	 * @param col
	 * @return A[row,:]*B[:,col]
	 */
	public double mult(SparseMatrixColMajor B, int row, int col) {
		Map<Integer,Double> rowA = m.get(row);
		Map<Integer,Double> colB = B.m.get(col);
		if(rowA==null || colB==null)
			return 0.0;
		double rlt = 0.0;
		if(rowA.size() < colB.size()) {
			for(Entry<Integer,Double> e : rowA.entrySet()) {
				int c = e.getKey();
				Double vB = colB.get(c);
				if(vB != null) rlt += e.getValue()*vB;
			}
		} else {
			for(Entry<Integer,Double> e : colB.entrySet()) {
				int r = e.getKey();
				Double vA = rowA.get(r);
				if(vA != null) rlt += e.getValue()*vA;
			}
		}
		return rlt;
	}
	
	/**
	 *Constructs and returns a new <tt>SparseVecotor</tt> view representing the columns of the given row.
	 *The returned view is backed by this matrix, so changes in the returned view are reflected in this matrix, and vice-versa.
	 * 
	 * @param row
	 */
	public SparseVectorHashMap viewRow(int row) {
		SparseVectorHashMap rlt = new SparseVectorHashMap(this.rowDim,
				this.m.get(row),false);
		return rlt;
	}
	
	/**
	 * Swap <tt>row1</tt> and <tt>row2</tt>
	 * 
	 * @param row1
	 * @param row2
	 */
	public void swapRow(int row1, int row2) {
		Map<Integer, Double> tmp = m.get(row1);
		m.put(row1, m.get(row2));
		m.put(row2, tmp);
	}
	
	/**
	 * Write this matrix to a file with Matlab mat file format.
	 * The variable name in matlab workspace is specified by <tt>setName()</tt>.
	 * Default variable name is <tt>"SparseMatrix"+UniqueSequenceNumber</tt>.
	 * <p>
	 * If more than one matrix need to be written in a single mat file use <tt>MatlabMatFileWriter</tt> instead.
	 * 
	 * @param fileName
	 */
	public void writeMatFile(String fileName) {
		MatlabMatFileWriter w = new MatlabMatFileWriter();
		w.addSparseMatrix(this);
		w.writeFile(fileName);
	}
	
	public void writeSimpleFile(String fileName) {
		throw new UnsupportedOperationException();
	}

	/**
	 * Return a row-major iterator
	 */
	@Override
	public Iterator<MatrixEntry> iterator() {
		return new SMIterator(this.m.entrySet().iterator());
	}
	
    /**
     * Iterator over this sparse matrix.
     */
    class SMIterator implements Iterator<MatrixEntry> {
        /**
         * Matrix cursor
         */
    	Iterator<Entry<Integer, Map<Integer,Double>>> rowIter;
    	Iterator<Entry<Integer, Double>> colIter;
   
    	SMIterator(Iterator<Entry<Integer, Map<Integer,Double>>> iter) {
    		rowIter = iter;
    		while(rowIter.hasNext()) {
	    		Entry<Integer, Map<Integer,Double>> nextRow = rowIter.next();
	    		colIter = nextRow.getValue().entrySet().iterator();
	    		if(colIter.hasNext()) {
	    			entry.row = nextRow;
	    			return;
	    		}
    		}
    	}
    	
        /**
         * Matrix entry
         */
        final SMEntry entry = new SMEntry();

        public boolean hasNext() {
    		if(colIter!=null && colIter.hasNext()) 
    			return true;
    		else {
        		while(rowIter.hasNext()) {
    	    		Entry<Integer, Map<Integer,Double>> nextRow = rowIter.next();
    	    		colIter = nextRow.getValue().entrySet().iterator();
    	    		if(colIter.hasNext()) {
    	    			entry.row = nextRow;
    	    			return true;
    	    		}
        		}
    		}
    		return false;
        }

        public MatrixEntry next() {
        	Entry<Integer, Double> ele = null;
        	if(colIter.hasNext()) {
        		ele = colIter.next();
        		entry.eleInRow = ele;
        	} else if(rowIter.hasNext()) {
        		Entry<Integer, Map<Integer,Double>> nextRow = rowIter.next();
        		entry.row = nextRow;
        		
        		colIter = nextRow.getValue().entrySet().iterator();
        		ele = colIter.next();
        		entry.eleInRow = ele;
        	}
            return entry;
        }

        public void remove() {
            colIter.remove();
        }

    }

    /**
     * Matrix entry backed by the matrix.
     */
    class SMEntry implements MatrixEntry {

    	private Entry<Integer, Map<Integer,Double>> row;
    	private Entry<Integer, Double> eleInRow;

		@Override
		public int getRow() {
			return row.getKey();
		}

		@Override
		public int getCol() {
			return eleInRow.getKey();
		}

		@Override
		public double getValue() {
			return eleInRow.getValue();
		}

		@Override
		public void setValue(double value) {
			eleInRow.setValue(value);
		}
    }

	@Override
	public double apply(int row, int col) {
		return this.get(row, col);
	}

	@Override
	public void update(int row, int col, double value) {
		this.set(row, col, value);
	}	
}
