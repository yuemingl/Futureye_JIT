package edu.uta.futureye.algebra;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Map.Entry;

import edu.uta.futureye.algebra.SparseVectorHashMap.SVEntry;
import edu.uta.futureye.algebra.intf.BlockVector;
import edu.uta.futureye.algebra.intf.SparseVector;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.algebra.intf.VectorEntry;
import edu.uta.futureye.io.MatlabMatFileWriter;
import edu.uta.futureye.util.FutureyeException;
import edu.uta.futureye.util.Sequence;

/**
 * Block vector
 * 
 * @author liuyueming
 *
 */
public class SparseBlockVector implements BlockVector<SparseVector>, SparseVector {
	protected int blockDim = 0;
	protected double defaultValue = 0.0;
	protected Map<Integer,SparseVector> data = 
		new HashMap<Integer,SparseVector>();
	protected String name = this.getClass().getSimpleName()+Sequence.getInstance().nextSeq();
	
	public SparseBlockVector(int blockDim) {
		this.blockDim = blockDim;
	}
	
	public SparseBlockVector(int blockDim, double defaultValue) {
		this.blockDim = blockDim;
		this.defaultValue = defaultValue;
	}

	@Override
	public void setDim(int dim) {
		throw new UnsupportedOperationException();
	}
	
	@Override
	public int getDim() {
		int rlt = 0;
		for(Entry<Integer,SparseVector> entry : data.entrySet()) 
			rlt += entry.getValue().getDim();
		return rlt;
	}
	
	@Override
	public void set(int index, double value) {
		int base = 0;
		Vector v = null;
		for(int bi=1;bi<=this.blockDim;bi++) {
			v = data.get(bi);
			int upper = v.getDim();
			if(base < index && index <= base + upper) {
				v.set(index - base, value);
				return;
			} else {
				base += upper;
			}
		}
	}
	
	@Override
	public SparseBlockVector set(Vector v) {
		if(v instanceof SparseBlockVector) {
			for(int bi=1;bi<=this.blockDim;bi++) {
				this.data.get(bi).set(((SparseBlockVector)v).getBlock(bi));
			}
		} else {
			for(int i=1;i<=v.getDim();i++) {
				this.set(i,v.get(i));
			}
		}
		return this;
	}
	
	@Override
	public SparseBlockVector add(double a, Vector v) {
		if(v instanceof SparseBlockVector) {
			for(int bi=1;bi<=this.blockDim;bi++) {
				this.data.get(bi).add(a,((SparseBlockVector)v).getBlock(bi));
			}
		} else {
			for(int i=1;i<=v.getDim();i++) {
				this.add(i,v.get(i));
			}
		}
		return this;
	}
	
	@Override
	public double get(int index) {
		int base = 0;
		Vector v = null;
		for(int bi=1;bi<=this.blockDim;bi++) {
			v = data.get(bi);
			int upper = v.getDim();
			if(base < index && index <= base + upper) {
				return v.get(index - base);
			} else {
				base += upper;
			}
		}
		throw new FutureyeException("index="+index);
	}

	@Override
	public void add(int index,double value) {
		set(index,get(index)+value);
	}
	
	@Override
	public SparseVector copy() {
		SparseBlockVector r = new SparseBlockVector(
				this.blockDim,this.defaultValue);
		for(Entry<Integer, SparseVector> e : data.entrySet()) {
			r.setBlock(e.getKey(), e.getValue().copy());
		}
		
		return r;
	}
	
	@Override
	public double norm2() {
		return Math.sqrt(this.dot(this));
	}
	
	@Override
	public double normInf() {
		Double max = Double.MIN_VALUE;
		Vector v = null;
		for(int bi=1;bi<=this.blockDim;bi++) {
			v = data.get(bi);
			double vInf = v.normInf();
			if(vInf > max) max = vInf;
		}
		return max;
	}
	
	@Override
	public double dot(Vector v2) {
		double rlt = 0.0;
		if(this.getDim() != v2.getDim()) {
			FutureyeException e = 
				new FutureyeException("ERROR: Vector dot product dim1="+this.getDim()+" != dim2="+v2.getDim());
			e.printStackTrace();
			System.exit(0);
		} else {
			for(int i=1;i<=getDim();i++) {
				rlt += this.get(i)*v2.get(i);
			}
		}
		return rlt;
	}

	@Override
	public int getBlockDim() {
		return this.blockDim;
	}

	@Override
	public Map<Integer, SparseVector> getAllBlock() {
		return this.data;
	}

	@Override
	public SparseVector getBlock(int index) {
		if(index > this.blockDim) 
			throw new FutureyeException(
					"index(="+index+") should <= blockDim(="+this.blockDim+")");
		return data.get(index);
	}

	@Override
	public void setBlock(int index, SparseVector v) {
		if(index > this.blockDim) 
			throw new FutureyeException(
					"index(="+index+") should <= blockDim(="+this.blockDim+")");
		data.put(index, v);
	}

	@Override
	public void print() {
		for(int i=1;i<=this.getDim();i++) {
			System.out.print(String.format("%8.4f", get(i))+"   ");
		}
		System.out.println();
	}

	@Override
	public void clearData() {
		for(Entry<Integer, SparseVector> e : this.data.entrySet()) {
			e.getValue().clearData();
		}
	}
	
	@Override
	public void clearAll() {
		for(Entry<Integer, SparseVector> e : this.data.entrySet()) {
			e.getValue().clearAll();
		}
		this.data.clear();
		this.blockDim = 0;
	}
	
	@Override
	public void setAll(double value) {
		for(Entry<Integer, SparseVector> e : this.data.entrySet()) {
			e.getValue().setAll(value);
		}
	}

	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append("SparseBlockVector(").
			append(this.blockDim).
			append("):N0R=").
			append(data.size()).
			append("\n");
		
		for(int i=1;i<=this.blockDim;i++) {
				sb.append("(").append(i).append(")=");
				sb.append(data.get(i)).append("\n");
		}
		return sb.toString();
	}

	@Override
	public SparseBlockVector add(Vector v) {
		throw new UnsupportedOperationException();
	}

	@Override
	public SparseBlockVector axpy(double a, Vector y) {
		throw new UnsupportedOperationException();
	}

	@Override
	public SparseBlockVector scale(double a) {
		throw new UnsupportedOperationException();
	}

	@Override
	public SparseBlockVector ax(double a) {
		throw new UnsupportedOperationException();
	}

	@Override
	public double norm1() {
		throw new UnsupportedOperationException();
	}

	@Override
	public SparseBlockVector set(double a, Vector v) {
		throw new UnsupportedOperationException();
	}

	@Override
	public SparseBlockVector shift(double dv) {
		throw new UnsupportedOperationException();
	}

	@Override
	public SparseBlockVector axDivy(double a, Vector y) {
		throw new UnsupportedOperationException();
	}

	@Override
	public SparseBlockVector axMuly(double a, Vector y) {
		throw new UnsupportedOperationException();
	}
	

	@Override
	public String getName() {
		return name;
	}

	@Override
	public SparseBlockVector setName(String name) {
		this.name = name;
		return this; 
	}

	@Override
	public int getNonZeroNumber() {
		int n=0;
		for(Entry<Integer,SparseVector> e : this.data.entrySet()) {
			n+=e.getValue().getNonZeroNumber();
		}
		return n;
	}

	@Override
	public Map<Integer, Double> getAll() {
		throw new UnsupportedOperationException();
	}

	@Override
	public SparseVector setAll(int nBase, Map<Integer, Double> dataMap) {
		throw new UnsupportedOperationException();
	}

	@Override
	public void writeMatFile(String fileName) {
		MatlabMatFileWriter w = new MatlabMatFileWriter();
		w.addSparseVector(this);
		w.writeFile(fileName);
	}

	@Override
	public void writeSimpleFile(String fileName) {
		throw new UnsupportedOperationException();
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
    	Iterator<VectorEntry> iter = null;
    	Iterator<Entry<Integer,SparseVector>> iterVector;
    	int nIndexBase = 0;
    	int nNextIndexBase = 0;
    	
        /**
         * Vector entry
         */
        final SVEntry entry = new SVEntry();
        
    	SVIterator(Iterator<Entry<Integer,SparseVector>> iter) {
    		this.iterVector = iter;
    		while(this.iterVector.hasNext()) {
    			Entry<Integer, SparseVector> v = this.iterVector.next();
    			nNextIndexBase = v.getValue().getDim();
    			this.iter = v.getValue().iterator();
    			if(this.iter.hasNext())
    				return;
    			nIndexBase += nNextIndexBase;
    		}
    	}
    	
        public boolean hasNext() {
        	if(iter !=null && iter.hasNext())
        		return true;
        	else {
        		while(this.iterVector.hasNext()) {
        			nIndexBase += nNextIndexBase;
        			Entry<Integer, SparseVector> v = this.iterVector.next();
        			nNextIndexBase = v.getValue().getDim();
        			this.iter = v.getValue().iterator();
        			if(this.iter.hasNext())
        				return true;
        		}
        	}
        	return false;
        }

        public VectorEntry next() {
        	entry.upate(nIndexBase, iter.next());
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

    	private int nIndexBase;
    	private VectorEntry e;

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
			return nIndexBase + e.getIndex();
		}
		
		private void upate(int nIndexBase, VectorEntry e) {
			this.nIndexBase = nIndexBase;
			this.e = e;
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
