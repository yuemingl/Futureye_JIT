package edu.uta.futureye.algebra;

import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.io.MatlabMatFileWriter;
import edu.uta.futureye.util.FutureyeException;
import edu.uta.futureye.util.Sequence;

/**
 * Space vector
 * <p>
 * 空间向量
 * 
 * @author liuyueming
 *
 */
public class SpaceVector implements Vector {
	protected int dim = 0;
	protected double[] data = null;
	protected String name = this.getClass().getSimpleName()+Sequence.getInstance().nextSeq();
	
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

	/**
	 * After setting dimension of the vector the old data will be lost
	 */
	@Override
	public void setDim(int dim) {
		this.dim = dim;
		data = new double[dim];
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
	public double get(int index) {
		return data[index-1];
	}
	
	@Override
	public void add(int index, double value) {
		set(index,get(index)+value);
	}
	
	/////////////////////////////////////////////////
	
	@Override
	public Vector set(Vector v) {
		for(int i=0;i<dim;i++)
			this.data[i] = v.get(i+1);
		return this;
	}

	@Override
	public Vector set(double a, Vector v) {
		for(int i=0;i<dim;i++)
			this.data[i] = a*v.get(i+1);
		return this;
	}
	
	@Override
	public Vector add(Vector v) {
		for(int i=0;i<dim;i++)
			this.data[i] += v.get(i+1);
		return this;
	}
	
	@Override
	public Vector add(double a, Vector v) {
		for(int i=0;i<dim;i++)
			this.data[i] += a*v.get(i+1);
		return this;
	}

	@Override
	public Vector scale(double a) {
		for(int i=0;i<dim;i++)
			this.data[i] = a*this.data[i];
		return this;
	}

	@Override
	public Vector ax(double a) {
		for(int i=0;i<dim;i++)
			this.data[i] = a*this.data[i];
		return this;
	}

	@Override
	public Vector axpy(double a, Vector y) {
		for(int i=0;i<dim;i++)
			this.data[i] = a*this.data[i] + y.get(i+1);
		return this;
	}
	
	@Override
	public Vector axMuly(double a, Vector y) {
		for(int i=0;i<dim;i++) {
			this.data[i] = a*this.data[i]*y.get(i+1);
		}
		return this;
	}
	
	@Override
	public Vector axDivy(double a, Vector y) {
		for(int i=0;i<dim;i++) {
			this.data[i] = a*this.data[i]/y.get(i+1);
		}
		return this;
	}
	
	@Override
	public Vector shift(double dv) {
		for(int i=0;i<dim;i++)
			this.data[i] = this.data[i]+dv;
		return this;
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
	public double norm1() {
		double rlt = 0.0;
		for(int i=0;i<dim;i++)
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
		for(int i=0;i<dim;i++) {
			double abs = Math.abs(data[i]);
			if(abs > max) max = abs;
		}
		return max;
	}

	/////////////////////////////////////////////////////
	
	@Override
	public Vector copy() {
		return new SpaceVector(this.data);
	}
	
	@Override
	public void print() {
		for(int i=1;i<=dim;i++) {
			System.out.print(get(i)+" ");
		}
		System.out.println("");
	}

	/////////////////////////////////////////////////
	/**
	 * cross product for 3D vectors a=(a1 a2 a3)' and b=(b1 b2 b3)'
	 * 叉乘（仅3维向量）
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
	///////////////////////////////////////////////////
	
	public String toString() {
		String rlt = "(";
		for(int i=0;i<dim-1;i++)
			rlt += data[i]+"  ";
		rlt += data[dim-1]+")";
		return rlt;
	}
	

	@Override
	public String getName() {
		return name;
	}

	@Override
	public Vector setName(String name) {
		this.name = name;
		return this; 
	}

	@Override
	public void setAll(double value) {
		for(int i=dim; --i>=0;)
			this.data[i] = value;
	}

	/**
	 * Write this vector to a file with Matlab mat file format.
	 * The variable name in matlab workspace is specified by <tt>setName()</tt>.
	 * Default variable name is <tt>"SpaceVector"+UniqueSequenceNumber</tt>.
	 * <p>
	 * If more than one vector need to be written in a single mat file use <tt>MatlabMatFileWriter</tt> instead.
	 * 
	 * @param fileName
	 */
	public void writeMatFile(String fileName) {
		MatlabMatFileWriter w = new MatlabMatFileWriter();
		w.addVector(this);
		w.writeFile(fileName);
	}

	@Override
	public void writeSimpleFile(String fileName) {
		// TODO Auto-generated method stub
		
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
