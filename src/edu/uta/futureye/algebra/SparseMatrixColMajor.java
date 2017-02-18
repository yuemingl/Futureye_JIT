/**
 * Copyright (c) 2010, nkliuyueming@gmail.com. All rights reserved.
 * 
 * 
 */
package edu.uta.futureye.algebra;

import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.MatrixEntry;
import edu.uta.futureye.algebra.intf.SparseMatrix;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.io.MatlabMatFileWriter;
import edu.uta.futureye.util.Sequence;

/**
 * Column-major map based storage sparse matrix
 * 
 * @author liuyueming
 *
 */
public class SparseMatrixColMajor implements SparseMatrix {
	protected int rowDim;
	protected int colDim;
	protected double defaultValue = 0.0;
	protected Map<Integer,Map<Integer,Double>> m = 
		new HashMap<Integer,Map<Integer,Double>>();
	protected String name = this.getClass().getSimpleName()+Sequence.getInstance().nextSeq();
	
	public SparseMatrixColMajor() {
	}
	
	public SparseMatrixColMajor(String name) {
		this.name = name;
	}	
	
	public SparseMatrixColMajor(int rowDim, int colDim) {
		this.rowDim = rowDim;
		this.colDim = colDim;
	}
	
	public SparseMatrixColMajor(String name, int rowDim, int colDim) {
		this.name = name;
		this.rowDim = rowDim;
		this.colDim = colDim;
	}
	
	public SparseMatrixColMajor(int rowDim, int colDim,double defaultValue) {
		this.rowDim = rowDim;
		this.colDim = colDim;
		this.defaultValue = defaultValue;
	}
	
	public SparseMatrixColMajor(String name,int rowDim, int colDim,double defaultValue) {
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
		Map<Integer,Double> aCol = m.get(col);
		if(aCol == null) {
			if(Math.abs(value) >= Matrix.zeroEps) {
				aCol = new HashMap<Integer,Double>();
				m.put(col, aCol);
				aCol.put(row, value);
			}
		} else {
			if(Math.abs(value) < Matrix.zeroEps)
				aCol.remove(row);
			else
				aCol.put(row, value);
		}
	}
	
	@Override
	public double get(int row, int col) {
		Map<Integer,Double> aCol = m.get(col);
		if(aCol == null) {
			return 0.0;
		} else {
			Double v = aCol.get(row);
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
			int nCol = rowEentry.getKey();
			Map<Integer, Double> col = rowEentry.getValue();
			for(Entry<Integer, Double> entry : col.entrySet()) {
				int nRow = entry.getKey();
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
		throw new UnsupportedOperationException();
	}
	
	@Override
	public SparseMatrix trans() {
		Map<Integer,Map<Integer,Double>> m2 = m;
		m = new HashMap<Integer,Map<Integer,Double>>();
		int dim = this.colDim;
		this.colDim = this.rowDim;
		this.rowDim = dim;
		for(Entry<Integer,Map<Integer,Double>> col : m2.entrySet()) {
			int nCol = col.getKey();
			for(Entry<Integer,Double> row : col.getValue().entrySet()) {
				int nRow = row.getKey();
				set(nCol, nRow, row.getValue());
			}
		}
		return this;
	}
	
	/**
	 * An overriding method can also return a subtype of the type returned by the overridden method.
	 * This is called a covariant return type.
	 */
	@Override
	public SparseMatrixColMajor copy() {
		SparseMatrixColMajor newM = new SparseMatrixColMajor(this.rowDim,this.colDim);
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
			"):N0C="+m.size();
	}

	////////////////////////////////////////////////////
	
	/**
	 * 返回列压缩存储方式的列索引数组，列号从0开始
	 */
	public int [][] getRowIndex() {
		int[][] rowIndex = new int[m.size()][];
		//for(int c=0; c<this.colDim; c++) {
		for(int c=this.colDim; --c>=0;) {
			Map<Integer,Double> col = m.get(c+1);
			rowIndex[c] = new int[col.size()];
			int r = 0;
			for(Entry<Integer,Double> row : col.entrySet()) {
				int nRow = row.getKey();
				rowIndex[c][r] = nRow-1;
				r++;
			}
		}
		return rowIndex;
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

	public int getNonZeroNumber() {
		int rlt = 0;
		for(Entry<Integer,Map<Integer,Double>> e1 : m.entrySet()) {
			rlt += e1.getValue().size();
		}
		return rlt;
	}
	
	/**
	 *Constructs and returns a new <tt>SparseVecotor</tt> view representing the row of the given column.
	 *The returned view is backed by this matrix, so changes in the returned view are reflected in this matrix, and vice-versa.
	 * 
	 * @param row
	 */
	public SparseVectorHashMap viewCol(int col) {
		SparseVectorHashMap rlt = new SparseVectorHashMap(this.rowDim,
				this.m.get(col),false);
		return rlt;
	}
	
	/**
	 * Return A[:,col]*vec, where A==this, vec.length=A.getRowDim()
	 * 
	 * @param vec
	 * @return
	 */
	public double multColumn(double[] vec, int col) {
		Map<Integer,Double> c = m.get(col);
		double rlt = 0.0;
		if(c != null) {
			for(Entry<Integer,Double> e : c.entrySet())
				rlt += e.getValue()*vec[e.getKey()-1];
		}
		return rlt;
	}
	
	/**
	 * Return A[:,col]*vec, where A==this, vec.length=A.getRowDim()
	 * 
	 * @param vec 
	 * @param vec Nonzero index list of vec, start form 1
	 * @return
	 */
	public double multColumn(double[] vec,List<Integer> nonzeroIndex, int col) {
		double rlt = 0.0;
		if(nonzeroIndex.size() == 0) 
			return rlt;
		Map<Integer,Double> c = m.get(col);
		if(c != null) {
			for(Integer idx : nonzeroIndex) {
				Double v = c.get(idx);
				if(v!=null)
					rlt += vec[idx-1]*v;
			}
		}
		return rlt;
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
    	Iterator<Entry<Integer, Map<Integer,Double>>> colIter;
    	Iterator<Entry<Integer, Double>> rowIter;
   
    	SMIterator(Iterator<Entry<Integer, Map<Integer,Double>>> iter) {
    		colIter = iter;
    		while(colIter.hasNext()) {
	    		Entry<Integer, Map<Integer,Double>> nextCol = colIter.next();
	    		rowIter = nextCol.getValue().entrySet().iterator();
	    		if(rowIter.hasNext()) {
	    			entry.col = nextCol;
	    			return;
	    		}
    		}
    	}
    	
        /**
         * Matrix entry
         */
        final SMEntry entry = new SMEntry();

        public boolean hasNext() {
    		if(rowIter!=null && rowIter.hasNext()) 
    			return true;
    		else {
        		while(colIter.hasNext()) {
    	    		Entry<Integer, Map<Integer,Double>> nextCol = colIter.next();
    	    		rowIter = nextCol.getValue().entrySet().iterator();
    	    		if(rowIter.hasNext()) {
    	    			entry.col = nextCol;
    	    			return true;
    	    		}
        		}
    		}
    		return false;
        }

        public MatrixEntry next() {
        	Entry<Integer, Double> ele = null;
        	if(rowIter.hasNext()) {
        		ele = rowIter.next();
        		entry.eleInCol = ele;
        	} else if(colIter.hasNext()) {
        		Entry<Integer, Map<Integer,Double>> nextCol = colIter.next();
        		entry.col = nextCol;
        		
        		rowIter = nextCol.getValue().entrySet().iterator();
        		ele = rowIter.next();
        		entry.eleInCol = ele;
        	}
            return entry;
        }

        public void remove() {
            rowIter.remove();
        }

    }

    /**
     * Matrix entry backed by the matrix.
     */
    class SMEntry implements MatrixEntry {

    	private Entry<Integer, Map<Integer,Double>> col;
    	private Entry<Integer, Double> eleInCol;

		@Override
		public int getRow() {
			return eleInCol.getKey();
		}

		@Override
		public int getCol() {
			return col.getKey();
		}

		@Override
		public double getValue() {
			return eleInCol.getValue();
		}

		@Override
		public void setValue(double value) {
			eleInCol.setValue(value);
		}
    }

	@Override
	public void writeSimpleFile(String fileName) {
		// TODO Auto-generated method stub
		
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
