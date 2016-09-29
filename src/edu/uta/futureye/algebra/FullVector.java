package edu.uta.futureye.algebra;

import edu.uta.futureye.algebra.intf.AlgebraVector;
import edu.uta.futureye.algebra.intf.Vector;

/**
 * Full(dense) vector. Data is stored in a double array.
 * 
 * @author liuyueming
 *
 */
public class FullVector implements AlgebraVector {
	/**
	 * Vector dimension(length)
	 */
	protected int dim = 0;
	
	/**
	 * Vector data stored in a double array.
	 */
	protected double[] data = null;

	public FullVector(int dim) {
		this.dim = dim;
		this.data = new double[this.dim];
		for(int i=dim; --i>=0;)
			data[i] = 0.0;
	}
	
	public FullVector(double[] data, boolean bCopy) {
		this.dim = data.length;
		if(bCopy) {
			this.data = new double[this.dim];
			System.arraycopy(data, 0, this.data, 0, this.dim);
		} else {
			this.data = data;
		}
	}
	
	public FullVector(int dim, double defaultValue) {
		this.dim = dim;
		this.data = new double[this.dim];
		for(int i=dim; --i>=0;)
			data[i] = defaultValue;
	}
	
	public FullVector(Vector v) {
		this.dim = v.getDim();
		this.data = new double[this.dim];
		for(int i=dim; --i>=0;)
			data[i] = v.get(i+1);
	}
	
	/**
	 * set every element to scale*Math.random()+shift
	 * 
	 * @param scale
	 * @param shift
	 */
	public void setRandom(double scale,double shift) {
		for(int i=dim; --i>=0;)
			data[i] = scale*Math.random()+shift;
	}
	
	@Override
	public int getDim() {
		return this.dim;
	}
	
	@Override
	public double[] getData() {
		return this.data;
	}
	
	@Override
	public FullVector set(AlgebraVector v) {
		double[] yData = v.getData();
		for(int i=dim; --i>=0;) {
			this.data[i] = yData[i];
		}
		return this;
	}

	@Override
	public FullVector set(double a, AlgebraVector v) {
		double[] yData = v.getData();
		for(int i=dim; --i>=0;) {
			this.data[i] = a*yData[i];
		}
		return this;
	}

	@Override
	public FullVector scale(double a) {
		for(int i=dim; --i>=0;)
			this.data[i] *= a;
		return this;
	}
	
	@Override
	public FullVector add(AlgebraVector v) {
		double[] yData = v.getData();
		for(int i=dim; --i>=0;) {
			this.data[i] += yData[i];
		}
		return this;
	}
	
	@Override
	public FullVector subtract(AlgebraVector v) {
		double[] yData = v.getData();
		for(int i=dim; --i>=0;) {
			this.data[i] -= yData[i];
		}
		return this;
	}
	
	@Override
	public FullVector add(double a, AlgebraVector v) {
		double[] yData = v.getData();
		for(int i=dim; --i>=0;) {
			this.data[i] = this.data[i] + a*yData[i];
		}
		return this;
	}

	@Override
	public double dot(AlgebraVector y) {
		double[] yData = y.getData();
		double rlt = 0.0;
		for(int i=dim; --i>=0;) {
			rlt += this.data[i]*yData[i];
		}
		return rlt;
	}

	@Override
	public FullVector ax(double a) {
		for(int i=dim; --i>=0;)
			this.data[i] *= a;
		return this;
	}

	@Override
	public FullVector axpy(double a, AlgebraVector y) {
		double[] yData = y.getData();
		for(int i=dim; --i>=0;) {
			this.data[i] = a*this.data[i] + yData[i];
		}
		return this;
	}
	
	@Override
	public FullVector axmy(double a, AlgebraVector y) {
		double[] yData = y.getData();
		for(int i=dim; --i>=0;) {
			this.data[i] = a*this.data[i] * yData[i];
		}
		return this;
	}

	@Override
	public double norm1() {
		double rlt = 0.0;
		for(int i=dim; --i>=0;)
			rlt += Math.abs(this.data[i]);
		return rlt;
	}	
	
	@Override
	public double norm2() {
		return Math.sqrt(this.dot(this));
	}
	
	@Override
	public double normInf() {
		Double max = Double.MIN_VALUE;
		for(int i=dim; --i>=0;) {
			double abs = Math.abs(data[i]);
			if(abs > max) max = abs;
		}
		return max;
	}

	@Override
	public void print() {
		for(int i=0;i<dim;i++) {
			System.out.print(String.format("%8.6f", this.data[i])+"   ");
		}
		System.out.println();
	}
	
	public class SparseData {
		public int[] index;
		public double[] data;
	}
	
	public SparseData getSparseData() {
		int n = 0;
		for(int i=0;i<this.dim;i++) {
			if(Double.compare(this.data[i], 0.0) != 0)
				n++;
		}
		SparseData sd = new SparseData();
		sd.index = new int[n];
		sd.data = new double[n];
		n=0;
		for(int i=0;i<this.dim;i++) {
			if(Double.compare(data[i], 0.0) != 0) {
				sd.index[n] = i;
				sd.data[n] = this.data[i];
				n++;
			}
		}
		return sd;
	}
	
	public SparseVectorHashMap getSparseVector() {
		SparseVectorHashMap rlt = new SparseVectorHashMap(this.dim);
		for(int i=dim; --i>=0;)
			rlt.set(i+1, this.data[i]);
		return rlt;
	}

	public FullVector copy() {
		int dim = this.dim;
		FullVector rlt = new FullVector(dim);
		for(int i=dim; --i>=0;) {
			rlt.data[i] = this.data[i];
		}
		return rlt;
	}
	
	/**
	 * 考虑该函数是否有必要
	 * @param index
	 * @return
	 */
	public double get(int index) {
		return data[index-1];
	}

}
