package edu.uta.futureye.function.basic;

import edu.uta.futureye.function.AbstractMathFun;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.util.Constant;

public class FXYZ extends AbstractMathFun{

	public FXYZ() {
		varNames.add(Constant.x);
		varNames.add(Constant.y);
		varNames.add(Constant.z);
	}
	
	public FXYZ(String varName) {
		super(varName);
	}

	
	@Override
	public MathFunc _d(String varName) {
		return null;
	}

	@Override
	public double apply(Variable v) {
		return v.get(getVarNames().get(0));
	}
	
	public String toString() {
		return getVarNames().get(0);
	}
}