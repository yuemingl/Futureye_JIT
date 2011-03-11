package edu.uta.futureye.function.basic;

import java.util.LinkedList;
import java.util.List;

import edu.uta.futureye.algebra.SpaceVector;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.function.intf.VectorFunction;
import edu.uta.futureye.function.operator.FOBasic;
import edu.uta.futureye.util.FutureyeException;
import edu.uta.futureye.util.Utils;

public class SpaceVectorFunction implements VectorFunction {
	protected int dim = 0;
	private List<String> varNames = null;
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
				//TODO copy ???
				data[i] = f[i];
				varNames = Utils.mergeList(varNames, f[i].varNames());
			}
		}
	}
	
	public SpaceVectorFunction(Vector v) {
		this.dim = v.getDim();
		data = new Function[dim];
		for(int i=0;i<dim;i++) {
			data[i] = new FC(v.get(i+1));
		}
	}

	@Override
	public VectorFunction copy() {
		VectorFunction rlt = new SpaceVectorFunction(this.data);
		//TODO rlt.setVarNames(varNames);
		return rlt;
	}

	@Override
	public Function dot(VectorFunction b) {
		if(dim != b.getDim()) {
			Exception e = new FutureyeException("Dims between two vector functions must be same!");
			e.printStackTrace();
			return null;
		}
		Function rlt = new FC(0);
		for(int i=0; i<dim; i++) {
			rlt = FOBasic.Plus(rlt, 
					FOBasic.Mult(data[i], b.get(i+1)));
		}
		return rlt;
	}
	
	/**
	 * 支持VectorFunction与Vector对象之间的内积
	 * @param b
	 * @return
	 */
	@Override
	public Function dot(Vector b) {
		return dot(new SpaceVectorFunction(b));
	}

	@Override
	public Function get(int index) {
		return data[index-1];
	}

	@Override
	public int getDim() {
		return dim;
	}

	@Override
	public Function norm2() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Function normInf() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void print() {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void set(int index, Function value) {
		data[index-1] = value; //copy???
		varNames = Utils.mergeList(varNames, value.varNames());
	}

	@Override
	public void setVarNames(List<String> varNames) {
		this.varNames = varNames;
	}

	@Override
	public Vector value(Variable v) {
		Vector rlt = new SpaceVector(dim);
		for(int i=0;i<dim;i++)
			rlt.set(i+1, data[i].value(v));
		return rlt;
	}

	@Override
	public List<String> varNames() {
		return varNames;
	}

	public String toString() {
		String rlt = "(";
		for(int i=0;i<dim;i++)
			rlt += data[i].toString()+"  ";
		return rlt+")";
	}
}
