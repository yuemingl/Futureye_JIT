package edu.uta.futureye.algebra;

import edu.uta.futureye.algebra.intf.AlgebraVector;

public class FullVector implements AlgebraVector {
	protected int dim = 0;
	protected double[] data = null;

	public FullVector(int dim) {
		this.dim = dim;
		this.data = new double[this.dim];
		for(int i=0;i<dim;i++)
			data[i] = 0.0;
	}
	
	public FullVector(int dim, double defaultValue) {
		this.dim = dim;
		this.data = new double[this.dim];
		for(int i=0;i<dim;i++)
			data[i] = defaultValue;
	}
	
	public FullVector(SparseVector v) {
		this.dim = v.getDim();
		this.data = new double[this.dim];
		for(int i=0;i<dim;i++)
			data[i] = v.get(i+1);
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
	public AlgebraVector set(AlgebraVector v) {
		double[] yData = v.getData();
		for(int i=0;i<this.dim;i++) {
			this.data[i] = yData[i];
		}
		return this;
	}
	
	@Override
	public AlgebraVector scale(double a) {
		for(int i=0;i<this.dim;i++)
			this.data[i] *= a;
		return this;
	}
	
	@Override
	public AlgebraVector plus(AlgebraVector v) {
		double[] yData = v.getData();
		for(int i=0;i<this.dim;i++) {
			this.data[i] += yData[i];
		}
		return this;
	}
	
	@Override
	public AlgebraVector minus(AlgebraVector v) {
		double[] yData = v.getData();
		for(int i=0;i<this.dim;i++) {
			this.data[i] -= yData[i];
		}
		return this;
	}
	
	@Override
	public AlgebraVector add(double a, AlgebraVector v) {
		double[] yData = v.getData();
		for(int i=0;i<this.dim;i++) {
			this.data[i] = this.data[i] + a*yData[i];
		}
		return this;
	}

	@Override
	public double dot(AlgebraVector y) {
		double[] yData = y.getData();
		double rlt = 0.0;
		for(int i=0; i<dim; i++) {
			rlt += this.data[i]*yData[i];
		}
		return rlt;
	}

	@Override
	public AlgebraVector ax(double a) {
		for(int i=0;i<this.dim;i++)
			this.data[i] *= a;
		return this;
	}

	@Override
	public AlgebraVector axpy(double a, AlgebraVector y) {
		double[] yData = y.getData();
		for(int i=0;i<this.dim;i++) {
			this.data[i] = a*this.data[i] + yData[i];
		}
		return this;
	}
	
	@Override
	public AlgebraVector axmy(double a, AlgebraVector y) {
		double[] yData = y.getData();
		for(int i=0;i<this.dim;i++) {
			this.data[i] = a*this.data[i] * yData[i];
		}
		return this;
	}

	@Override
	public double norm2() {
		return Math.sqrt(this.dot(this));
	}
	
	@Override
	public double normInf() {
		Double max = Double.MIN_VALUE;
		for(double d : data) {
			if(d > max) max = d;
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

}
