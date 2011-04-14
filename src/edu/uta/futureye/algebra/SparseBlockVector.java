package edu.uta.futureye.algebra;

import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;

import edu.uta.futureye.algebra.intf.BlockVector;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.util.FutureyeException;

public class SparseBlockVector implements BlockVector {
	protected int blockDim = 0;
	protected double defaultValue = 0.0;
	protected Map<Integer,Vector> data = 
		new HashMap<Integer,Vector>();
	
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
		for(Entry<Integer,Vector> entry : data.entrySet()) 
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
	public Vector set(Vector v) {
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
	public Vector add(double a, Vector v) {
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
	public Vector copy() {
		SparseBlockVector r = new SparseBlockVector(
				this.blockDim,this.defaultValue);
		for(Entry<Integer, Vector> e : data.entrySet()) {
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
	public Map<Integer, Vector> getAllBlock() {
		return this.data;
	}

	@Override
	public Vector getBlock(int index) {
		return data.get(index);
	}

	@Override
	public void setBlock(int index, Vector v) {
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
	public void clear() {
		throw new UnsupportedOperationException();
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
	public Vector add(Vector v) {
		throw new UnsupportedOperationException();
	}

	@Override
	public Vector axpy(double a, Vector y) {
		throw new UnsupportedOperationException();
	}

	@Override
	public Vector scale(double a) {
		throw new UnsupportedOperationException();
	}

	@Override
	public Vector ax(double a) {
		throw new UnsupportedOperationException();
	}

	@Override
	public double norm1() {
		throw new UnsupportedOperationException();
	}

	@Override
	public Vector set(double a, Vector v) {
		throw new UnsupportedOperationException();
	}

	@Override
	public Vector shift(double dv) {
		throw new UnsupportedOperationException();
	}
}
