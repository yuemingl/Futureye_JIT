package edu.uta.futureye.function.basic;

import java.util.LinkedList;

import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.function.AbstractVectorFunction;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.function.intf.VectorFunction;
import edu.uta.futureye.util.FutureyeException;
import edu.uta.futureye.util.Utils;

public class SpaceVectorFunction extends AbstractVectorFunction {
	protected Function[] data = null;
	
	public SpaceVectorFunction(int dim) {
		this.dim = dim;
		data = new Function[dim];
	}
	
	public SpaceVectorFunction(Function ...f) {
		if(f == null || f.length ==0) {
			Exception e = new FutureyeException("Dim of SpaceVectorFunction should be > 0!");
			e.printStackTrace();
			return;
		} else {
			dim = f.length;
			data = new Function[dim];
			varNames = new LinkedList<String>();
			for(int i=0; i<f.length; i++) {
				data[i] = f[i];
				varNames = Utils.mergeList(varNames, f[i].varNames());
			}
		}
	}
	
	public SpaceVectorFunction(Vector v) {
		this.dim = v.getDim();
		data = new Function[dim];
		for(int i=1; i<=dim; i++) {
			data[i-1] = new FC(v.get(i));
		}
	}
	
	@Override
	public Function get(int index) {
		return data[index-1];
	}

	@Override
	public void set(int index, Function value) {
		data[index-1] = value;
		varNames = Utils.mergeList(varNames, value.varNames());
	}
	
	/////////////////////////////////////////////////////
	
	@Override
	public VectorFunction set(VectorFunction v) {
		for(int i=1; i<=dim; i++) {
			data[i-1] = v.get(i);
		}
		this.setVarNames(v.varNames());
		return this;
	}
	
	@Override
	public VectorFunction set(double a, VectorFunction v) {
		for(int i=1; i<=dim; i++) {
			data[i-1] = v.get(i).M(a);
		}
		this.setVarNames(v.varNames());
		return this;
	}
	
	@Override
	public VectorFunction add(VectorFunction v) {
		for(int i=1; i<=dim; i++) {
			data[i-1] = data[i-1].A(v.get(i));
		}
		varNames = Utils.mergeList(varNames, v.varNames());
		return this;
	}
	
	@Override
	public VectorFunction add(double a, VectorFunction g) {
		for(int i=1; i<=dim; i++) {
			data[i-1] = data[i-1].A(g.get(i).M(a));
		}
		varNames = Utils.mergeList(varNames, g.varNames());
		return this;
	}
	
	@Override
	public VectorFunction scale(double a) {
		for(int i=1; i<=dim; i++) {
			data[i-1] = data[i-1].M(a);
		}
		return this;
	}
	
	@Override
	public VectorFunction ax(double a) {
		for(int i=1; i<=dim; i++) {
			data[i-1] = data[i-1].M(a);
		}
		return this;
	}
	
	@Override
	public VectorFunction axpy(double a, VectorFunction g) {
		for(int i=1; i<=dim; i++) {
			data[i-1] = data[i-1].M(a).A(g.get(i));
		}
		varNames = Utils.mergeList(varNames, g.varNames());
		return this;
	}
	
	/////////////////////////////////////////////////////
	
	@Override
	public VectorFunction copy() {
		SpaceVectorFunction rlt = new SpaceVectorFunction(dim);
		varNames = new LinkedList<String>();
		for(int i=0; i<=dim; i++) {
			rlt.data[i] = this.data[i].copy();
			varNames = Utils.mergeList(varNames, this.data[i].varNames());
		}
		return rlt;
	}

	@Override
	public void print() {
		System.out.println(this.toString());
	}
	
	@Override
	public String toString() {
		String rlt = "(";
		for(int i=0;i<dim;i++)
			rlt += data[i].toString()+"  ";
		return rlt+")";
	}
}
