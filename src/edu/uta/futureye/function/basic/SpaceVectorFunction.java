package edu.uta.futureye.function.basic;

import java.util.LinkedList;

import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.function.VectorMathFuncBase;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.function.intf.VectorMathFunc;
import edu.uta.futureye.util.FutureyeException;
import edu.uta.futureye.util.Utils;

public class SpaceVectorFunction extends VectorMathFuncBase {
	protected MathFunc[] data = null;
	
	public SpaceVectorFunction(int dim) {
		this.dim = dim;
		data = new MathFunc[dim];
	}
	
	public SpaceVectorFunction(MathFunc ...f) {
		if(f == null || f.length ==0) {
			Exception e = new FutureyeException("Dim of SpaceVectorFunction should be > 0!");
			e.printStackTrace();
			return;
		} else {
			dim = f.length;
			data = new MathFunc[dim];
			varNames = new LinkedList<String>();
			for(int i=0; i<f.length; i++) {
				data[i] = f[i];
				varNames = Utils.mergeList(varNames, f[i].getVarNames());
			}
		}
	}
	
	public SpaceVectorFunction(Vector v) {
		this.dim = v.getDim();
		data = new MathFunc[dim];
		for(int i=1; i<=dim; i++) {
			data[i-1] = new FC(v.get(i));
		}
	}
	
	@Override
	public MathFunc get(int index) {
		return data[index-1];
	}

	@Override
	public void set(int index, MathFunc value) {
		data[index-1] = value;
		varNames = Utils.mergeList(varNames, value.getVarNames());
	}
	
	/////////////////////////////////////////////////////
	
	@Override
	public VectorMathFunc set(VectorMathFunc v) {
		for(int i=1; i<=dim; i++) {
			data[i-1] = v.get(i);
		}
		this.setVarNames(v.varNames());
		return this;
	}
	
	@Override
	public VectorMathFunc set(double a, VectorMathFunc v) {
		for(int i=1; i<=dim; i++) {
			data[i-1] = v.get(i).M(a);
		}
		this.setVarNames(v.varNames());
		return this;
	}
	
	@Override
	public VectorMathFunc add(VectorMathFunc v) {
		for(int i=1; i<=dim; i++) {
			data[i-1] = data[i-1].A(v.get(i));
		}
		varNames = Utils.mergeList(varNames, v.varNames());
		return this;
	}
	
	@Override
	public VectorMathFunc add(double a, VectorMathFunc g) {
		for(int i=1; i<=dim; i++) {
			data[i-1] = data[i-1].A(g.get(i).M(a));
		}
		varNames = Utils.mergeList(varNames, g.varNames());
		return this;
	}
	
	@Override
	public VectorMathFunc scale(double a) {
		for(int i=1; i<=dim; i++) {
			data[i-1] = data[i-1].M(a);
		}
		return this;
	}
	
	@Override
	public VectorMathFunc ax(double a) {
		for(int i=1; i<=dim; i++) {
			data[i-1] = data[i-1].M(a);
		}
		return this;
	}
	
	@Override
	public VectorMathFunc axpy(double a, VectorMathFunc g) {
		for(int i=1; i<=dim; i++) {
			data[i-1] = data[i-1].M(a).A(g.get(i));
		}
		varNames = Utils.mergeList(varNames, g.varNames());
		return this;
	}
	
	/////////////////////////////////////////////////////
	
	@Override
	public VectorMathFunc copy() {
		SpaceVectorFunction rlt = new SpaceVectorFunction(dim);
		varNames = new LinkedList<String>();
		for(int i=0; i<=dim; i++) {
			rlt.data[i] = this.data[i].copy();
			varNames = Utils.mergeList(varNames, this.data[i].getVarNames());
		}
		return rlt;
	}

	@Override
	public String getExpression() {
		String rlt = "[";
		for(int i=0;i<dim;i++) {
			rlt += data[i].toString();
			if(i<dim-1)
				rlt += ", ";
		}
		return rlt+"]";
	}
}
