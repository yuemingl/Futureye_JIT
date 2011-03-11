package edu.uta.futureye.function.basic;

import edu.uta.futureye.function.AbstractFunction;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.intf.Function;

public class FXYZ extends AbstractFunction{

	public FXYZ() {
	}
	
	public FXYZ(String varName) {
		super(varName);
	}

	
	@Override
	public Function d(String varName) {
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