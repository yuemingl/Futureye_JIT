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

import edu.uta.futureye.algebra.intf.SparseVector;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.algebra.intf.VectorEntry;
import edu.uta.futureye.io.MatlabMatFileWriter;
import edu.uta.futureye.util.FutureyeException;
import edu.uta.futureye.util.Sequence;

/**
 * HashMap based storage sparse vector implementation
 * <p>
 * This implementation provides constant-time performance for the basic operations (get and set), 
 * assuming the background hash function disperses the elements properly among the buckets.
 * 
 * @author liuyueming
 *
 */
public class SparseVectorHashMap implements SparseVector {
	protected int dim = 0;
	protected double defaultValue = 0.0;
	protected Map<Integer,Double> data = 
		new HashMap<Integer,Double>();
	protected String name = this.getClass().getSimpleName()+Sequence.getInstance().nextSeq();

	public SparseVectorHashMap() {
	}
	
	public SparseVectorHashMap(int dim) {
		this.dim = dim;
	}
	
	public SparseVectorHashMap(int dim, double defaultValue) {
		this.dim = dim;
		this.defaultValue = defaultValue;
	}
	
	/**
	 * Constructs a sparse vector with the given parameters with zero tolerance.
	 * The dimension of the vector is the same as the number of parameters
	 * 
	 * @param a
	 */
	public SparseVectorHashMap(double ...a) {
		if(a == null || a.length ==0) {
			throw new FutureyeException("Dim of SparseVector should be > 0!");
		} else {
			dim = a.length;
			for(int i=0; i<a.length; i++) {
				if(Double.compare(a[i], 0.0) != 0)
					data.put(i+1, a[i]);
			}
		}
	}
	
	/**
	 * Constructs a sparse vector with the given map <tt>data</tt>. 
	 * If <tt>bCopy==false</tt> the receiver is backed by map <tt>data</tt>, 
	 * so changes in the receiver are reflected in the map <tt>data</tt>, and vice-versa.
	 * 
	 * @param dim
	 * @param data
	 * @param bCopy
	 */
	public SparseVectorHashMap(int dim, Map<Integer,Double> data, boolean bCopy) {
		this.dim = dim;
		if(bCopy) {
			this.data.putAll(data);
		} else {
			this.data = data;
		}
	}
	
	/**
	 * Reset dimension will cause data loss!
	 */
	@Override
	public void setDim(int dim) {
		this.dim = dim;
		this.clearData();
	}
	
	@Override
	public int getDim() {
		return dim;
	}

	@Override
	public void set(int index, double value) {
		if(index>dim) {
			throw new FutureyeException("index("+index+") > dim("+dim+")");
		}
		data.put(index, value);
	}
	
	@Override
	public SparseVectorHashMap set(Vector v) {
		this.dim = v.getDim();
		if(v instanceof SparseVectorHashMap) {
			SparseVectorHashMap tmp = (SparseVectorHashMap)v;
			this.data.clear();
			for(Entry<Integer, Double> e : tmp.data.entrySet()) {
				this.data.put(e.getKey(), e.getValue());
			}
			this.defaultValue = tmp.defaultValue;
		} else {
			for(int i=1;i<=v.getDim();i++) {
				this.set(i,v.get(i));
			}
		}
		return this;
	}
	
	@Override
	public SparseVectorHashMap set(double a, Vector v) {
		this.dim = v.getDim();
		if(v instanceof SparseVectorHashMap) {
			SparseVectorHashMap tmp = (SparseVectorHashMap)v;
			this.data.clear();
			for(Entry<Integer, Double> e : tmp.data.entrySet()) {
				this.data.put(e.getKey(), a*e.getValue());
			}
			this.defaultValue = a*tmp.defaultValue;
		} else {
			for(int i=1;i<=v.getDim();i++) {
				this.set(i,a*v.get(i));
			}
		}
		return this;
	}
	
	@Override
	public double get(int index) {
		if(index>dim) {
			throw new FutureyeException("index("+index+") > dim("+dim+")");
		}
		Double val = data.get(index);
		if(val == null) {
			return this.defaultValue;
		} else {
			return val;
		}
	}

	@Override
	public void add(int index,double value) {
		set(index,get(index)+value);
	}
	
	@Override
	public SparseVectorHashMap add(Vector v) {
		if(v instanceof SparseVectorHashMap) {
			SparseVectorHashMap tmp = (SparseVectorHashMap)v;
			for(Entry<Integer, Double> e : tmp.data.entrySet()) {
				this.add(e.getKey(), e.getValue());
			}
			this.defaultValue += tmp.defaultValue;

		} else {
			for(int i=1;i<=v.getDim();i++) {
				this.add(i,v.get(i));
			}
		}
		return this;
	}
	
	@Override
	public SparseVectorHashMap add(double a, Vector v) {
		if(v instanceof SparseVectorHashMap) {
			SparseVectorHashMap tmp = (SparseVectorHashMap)v;
			for(Entry<Integer, Double> e : tmp.data.entrySet()) {
				this.add(e.getKey(), a*e.getValue());
			}
			this.defaultValue += a*tmp.defaultValue;

		} else {
			for(int i=1;i<=v.getDim();i++) {
				this.add(i,a*v.get(i));
			}
		}
		return this;
	}
	
	@Override
	public SparseVectorHashMap ax(double a) {
		for(Entry<Integer, Double> e : data.entrySet()) {
			this.set(e.getKey(), a*e.getValue());
		}
		this.defaultValue *= a;
		return this;
	}
	
	@Override
	public SparseVectorHashMap axpy(double a, Vector y) {
		this.scale(a).add(y);
		return this;
	}

	@Override
	public SparseVectorHashMap axMuly(double a, Vector y) {
		for(Entry<Integer, Double> e : data.entrySet()) {
			this.set(e.getKey(), a*e.getValue()*y.get(e.getKey()));
		}
		return this;
	}
	
	@Override
	public SparseVectorHashMap axDivy(double a, Vector y) {
		for(Entry<Integer, Double> e : data.entrySet()) {
			this.set(e.getKey(), a*e.getValue()/y.get(e.getKey()));
		}
		return this;
	}
	
	@Override
	public SparseVectorHashMap scale(double a) {
		for(Entry<Integer, Double> e : data.entrySet()) {
			this.set(e.getKey(), a*e.getValue());
		}
		this.defaultValue *= a;
		return this;
	}

	@Override
	public SparseVectorHashMap shift(double dv) {
		for(Entry<Integer, Double> e : data.entrySet()) {
			this.set(e.getKey(), e.getValue() + dv);
		}
		this.defaultValue += dv;
		return this;
	}

	@Override
	public double norm1() {
		throw new UnsupportedOperationException();
	}

	@Override
	public double norm2() {
		return Math.sqrt(this.dot(this));
	}
	
	@Override
	public double normInf() {
		Double max = Double.MIN_VALUE;
		for(Double d : data.values()) {
			double abs = Math.abs(d);
			if(abs > max) max = abs;
		}
		return max;
	}
	
	@Override
	public double dot(Vector v2) {
		double rlt = 0.0;
		if(this.getDim() != v2.getDim()) {
			throw new FutureyeException(
					"ERROR: Vector dot product dim1="+
					this.getDim()+" != dim2="+v2.getDim());
		} else {
			for(int i=1;i<=getDim();i++) {
				rlt += this.get(i)*v2.get(i);
			}
		}
		return rlt;
	}

	/**
	 * An overriding method can also return a subtype of the type returned by the overridden method. 
	 * This is called a covariant return type.
	 */
	@Override
	public SparseVectorHashMap copy() {
		SparseVectorHashMap r = new SparseVectorHashMap(this.dim,this.defaultValue);
		for(Entry<Integer, Double> e : data.entrySet()) {
			r.set(e.getKey(), e.getValue());
		}
		return r;
	}
	
	@Override
	public void clearData() {
		this.data.clear();
	}
	
	@Override
	public void clearAll() {
		this.dim = 0;
		this.data.clear();
	}
	
	@Override
	public Vector setAll(double value) {
		for(Entry<Integer, Double> e : this.data.entrySet()) {
			e.setValue(value);
		}
		return this;
	}
	
	@Override
	public void print() {
		for(int i=1;i<=dim;i++) {
			System.out.print(String.format("%8.6f", get(i))+"   ");
		}
		System.out.println();
		System.out.println();
	}
	
	public String toString() {
		return "SparseVector("+
			this.dim+
			"):N0R="+data.size();
	}
	
	/////////////////////////////////////////////////

	/**
	 * Get all non zero values in this vector
	 * 
	 * @return Integer-double map containing non zero indices and values in this vector
	 */
	public Map<Integer,Double> getAll() {
		return this.data;
	}

	/**
	 * Get number of non zero values
	 * 
	 */
	public int getNonZeroNumber() {
		return this.data.size();
	}

	@Override
	public String getName() {
		return name;
	}

	@Override
	public SparseVectorHashMap setName(String name) {
		this.name = name;
		return this; 
	}
	
	/**
	 * Write this vector to a file with Matlab mat file format.
	 * The variable name in matlab workspace is specified by <tt>setName()</tt>.
	 * Default variable name is <tt>"SparseVector"+UniqueSequenceNumber</tt>.
	 * <p>
	 * If more than one vector need to be written in a single mat file use <tt>MatlabMatFileWriter</tt> instead.
	 * 
	 * @param fileName
	 */
	public void writeMatFile(String fileName) {
		MatlabMatFileWriter w = new MatlabMatFileWriter();
		w.addSparseVector(this);
		w.writeFile(fileName);
	}

	@Override
	public void writeSimpleFile(String fileName) {
		// TODO Auto-generated method stub
		
	}
	
	@Override
	public SparseVector setAll(int nBase, Map<Integer, Double> dataMap) {
		for(Entry<Integer, Double> e : dataMap.entrySet()) {
			this.data.put(nBase+e.getKey(), e.getValue());
		}
		return this;
	}
	
	@Override
	public Iterator<VectorEntry> iterator() {
		return new SVIterator(this.data.entrySet().iterator());
	}
	
    /**
     * Iterator over this sparse vector.
     */
    class SVIterator implements Iterator<VectorEntry> {
       /**
         * Vector cursor
         */
    	Iterator<Entry<Integer,Double>> iter;
   
        /**
         * Vector entry
         */
        final SVEntry entry = new SVEntry();
        
     	SVIterator(Iterator<Entry<Integer,Double>> iter) {
    		this.iter = iter;
    	}
        public boolean hasNext() {
            return iter.hasNext();
        }

        public SVEntry next() {
        	entry.e = iter.next();
            return entry;
        }

        public void remove() {
        	iter.remove();
        }

    }

    /**
     * Vector entry backed by the vector.
     */
    class SVEntry implements VectorEntry {

    	private Entry<Integer, Double> e;

		@Override
		public double getValue() {
			return e.getValue();
		}

		@Override
		public void setValue(double value) {
			e.setValue(value);
		}

		@Override
		public int getIndex() {
			return e.getKey();
		}
    }

	@Override
	public double apply(int index) {
		return this.get(index);
	}

	@Override
	public void update(int index, double value) {
		this.set(index,value);
	}
}
