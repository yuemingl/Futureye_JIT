package edu.uta.futureye.algebra;

import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.util.FutureyeException;

public class SpaceVector implements Vector {
	protected int dim = 0;
	protected double[] data = null;
	
	public SpaceVector() {
	}
	
	public SpaceVector(int dim) {
		this.dim = dim;
		data = new double[dim];
	}
	
	public SpaceVector(double ...a) {
		if(a == null || a.length ==0) {
			Exception e = new FutureyeException("Dim of SpaceVector should be > 0!");
			e.printStackTrace();
			System.exit(0);
		} else {
			dim = a.length;
			data = new double[dim];
			for(int i=0; i<a.length; i++)
				data[i] = a[i];
		}
	}

	@Override
	public void setDim(int dim) {
		this.dim = dim;
	}
	
	@Override
	public int getDim() {
		return dim;
	}
	
	@Override
	public void set(int index, double value) {
		data[index-1] = value;
	}
	
	@Override
	public void set(Vector v) {
		SpaceVector tmp = (SpaceVector)v;
		for(int i=0;i<dim;i++)
			this.data[i] = tmp.data[i];
	}

	@Override
	public double get(int index) {
		return data[index-1];
	}
	
	@Override
	public void add(int index, double value) {
		set(index,get(index)+value);
	}
	
	@Override
	public void add(double a, Vector v) {
		SpaceVector tmp = (SpaceVector)v;
		for(int i=0;i<dim;i++)
			this.data[i] += a*tmp.data[i];
	}

	
	@Override
	public Vector copy() {
		return new SpaceVector(this.data);
	}

	@Override
	public double dot(Vector u) {
		double rlt = 0.0;
		if(dim != u.getDim()) {
			Exception e = new FutureyeException("Dims between two vectors must be same!");
			e.printStackTrace();
			System.exit(0);
		} else {
			for(int i=0;i<dim;i++) {
				rlt += data[i]*u.get(i+1);
			}
		}
		return rlt;
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
	public void clear() {
		this.dim = 0;
		this.data = null;
	}
	
	@Override
	public void print() {
		for(int i=1;i<=dim;i++) {
			System.out.print(String.format("%8.6f", get(i))+"   ");
		}
		System.out.println("");
	}

	/////////////////////////////////////////////////
	
	/**
	 * ²æ³Ë£¨½ö3Î¬ÏòÁ¿£©
	 * cross product for 3D vectors a=(a1 a2 a3)' and b=(b1 b2 b3)'
	 * 
	 *     |i   j  k|
	 * 3D  |a1 a2 a3|
	 *     |b1 b2 b3|
	 * 
	 *     |a2 a3|     |a1 a3|     |a1 a2|
	 *   = |b2 b3|*i - |b1 b3|*j + |b1 b2|*k
	 * 
	 * @param a
	 * @param b
	 * @return
	 */
	public SpaceVector crossProduct(SpaceVector b){
		SpaceVector r = new SpaceVector(3);
		if(dim!=3 || b.getDim()!=3 || dim != b.getDim()) {
			Exception e = new FutureyeException("Cross produce: all dims must be 3!");
			e.printStackTrace();
			return null;
		} else {
			r.set(1, get(2)*b.get(3)-get(3)*b.get(2));
			r.set(2, -get(1)*b.get(3)+get(3)*b.get(1));
			r.set(3, get(1)*b.get(2)-get(2)*b.get(1));
		}
		return r;
	}

	public static Vector ax(double a, Vector x) {
		int dim = x.getDim();
		Vector rlt = new SpaceVector(dim);
		for(int i=1;i<=dim;i++) {
			rlt.set(i, a*x.get(i));
		}
		return rlt;	
	}
	
	public static Vector axpy(double a, Vector x, Vector y) {
		int dim = x.getDim();
		Vector rlt = new SpaceVector(dim);
		for(int i=1;i<=dim;i++) {
			rlt.set(i, a*x.get(i)+y.get(i));
		}
		return rlt;
	}
	
	public static Vector axmy(double a, Vector x, Vector y) {
		int dim = x.getDim();
		Vector rlt = new SpaceVector(dim);
		for(int i=1;i<=dim;i++) {
			rlt.set(i, a*x.get(i)*y.get(i));
		}
		return rlt;
	}
	
	public String toString() {
		String rlt = "(";
		for(int i=0;i<dim;i++)
			rlt += data[i]+"  ";
		return rlt+")";
	}
}
