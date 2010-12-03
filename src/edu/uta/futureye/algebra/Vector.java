package edu.uta.futureye.algebra;

import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;

public class Vector {
	protected int dim;
	protected Map<Integer,Double> v = new HashMap<Integer,Double>();
	
	public Vector(int dim) {
		this.dim = dim;
	}
	
	public Vector(int dim, double initValue) {
		this.dim = dim;
		for(int i=1;i<=dim;i++) {
			set(i,initValue);
		}
	}
	
	public int getDim() {
		return dim;
	}
	
	public void set(int index, double value) {
		v.put(index, value);
	}
	
	public double get(int index) {
		Double val = v.get(index);
		if(val == null) {
			return 0.0;
		} else {
			return val;
		}
	}

	public void plusValue(int index,double value) {
		set(index,get(index)+value);
	}
	
	public Vector copy() {
		Vector r = new Vector(this.dim);
		for(Entry<Integer, Double> e : v.entrySet()) {
			r.set(e.getKey(), e.getValue());
		}
		return r;
	}
	
	public double norm2() {
		return Math.sqrt(this.dot(this));
	}
	
	public double normInf() {
		Double max = Double.MIN_VALUE;
		for(Double d : v.values()) {
			if(d > max) max = d;
		}
		return max;
	}
	
	public double dot(Vector v2) {
		double rlt = 0.0;
		if(this.getDim() != v2.getDim()) {
			System.out.println("ERROR: Vector dot product dim1="+this.getDim()+" != dim2="+v2.getDim());
		} else {
			for(int i=1;i<=getDim();i++) {
				rlt += this.get(i)*v2.get(i);
			}
		}
		return rlt;
	}
	
	public static double dot(Vector v1, Vector v2) {
		return v1.dot(v2);
	}
	
	/**
	 * 叉乘
	 * |i   j  k|
	 * |a1 a2 a3|
	 * |b1 b2 b3|
	 * 
	 *   |a2 a3|     |a1 a3|     |a1 a2|
	 * = |b2 b3|*i - |b1 b3|*j + |b1 b2|*k
	 * 
	 * @param a
	 * @param b
	 * @return
	 */
	public static Vector crossProduct3D(Vector a, Vector b) {
		Vector r = new Vector(3);
		if(a.getDim() != b.getDim()) {
			System.out.println("ERROR: Vector cross product dim1="+a.getDim()+" != dim2="+b.getDim());
		} else {
			r.set(1, a.get(2)*b.get(3)-a.get(3)*b.get(2));
			r.set(2, -a.get(1)*b.get(3)+a.get(3)*b.get(1));
			r.set(3, a.get(1)*b.get(2)-a.get(2)*b.get(1));
		}
		return r;
	}
	
	public Vector crossProduct3D(Vector b) {
		return crossProduct3D(this,b);
	}
	
	/**
	 * 快速创建一个Vector
	 * @param a1
	 * @param an
	 * @return
	 */
	public static Vector createVector(double a1, double ...an) {
		Vector r = null;
		if(an == null) {
			r = new Vector(1);
			r.set(1, a1);
		} else {
			r = new Vector(1+an.length);
			r.set(1, a1);
			for(int i=0;i<an.length;i++) {
				r.set(i+2, an[i]);
			}
		}
		return r;
	}
	
	public static Vector ax(double a, Vector x) {
		Vector rlt = new Vector(x.dim);
		for(int i=1;i<=x.dim;i++) {
			rlt.set(i, a*x.get(i));
		}
		return rlt;	
	}
	
	public static Vector axpy(double a, Vector x, Vector y) {
		Vector rlt = new Vector(x.dim);
		for(int i=1;i<=x.dim;i++) {
			rlt.set(i, a*x.get(i)+y.get(i));
		}
		return rlt;
	}
	
	public static Vector axmy(double a, Vector x, Vector y) {
		Vector rlt = new Vector(x.dim);
		for(int i=1;i<=x.dim;i++) {
			rlt.set(i, a*x.get(i)*y.get(i));
		}
		return rlt;
	}
	
	public void print() {
		for(int i=1;i<=dim;i++) {
			System.out.print(String.format("%8.6f", get(i))+"   ");
		}
		System.out.println("");
	}	
}
