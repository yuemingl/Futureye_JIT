package edu.uta.futureye.function.basic;

import edu.uta.futureye.function.AbstractFunction;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.util.Constant;

public class FXYZ extends AbstractFunction{

	public FXYZ() {
		varNames.add(Constant.x);
		varNames.add(Constant.y);
		varNames.add(Constant.z);
	}
	
	public FXYZ(String varName) {
		super(varName);
	}

	
	@Override
	public Function _d(String varName) {
		return null;
	}

	@Override
	public double value(Variable v) {
		return v.get(varNames().get(0));
	}
	
	public String toString() {
		return varNames().get(0);
	}
}