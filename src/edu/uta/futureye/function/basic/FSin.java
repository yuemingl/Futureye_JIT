package edu.uta.futureye.function.basic;

import java.util.List;

import edu.uta.futureye.function.AbstractMathFun;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.intf.MathFun;

public class FSin extends AbstractMathFun {
	
	public FSin() {
		super("x");
	}
	
	public FSin(String varName) {
		super(varName);
	}
	
	public FSin(List<String> varNames) {
		super(varNames);
	}
	
	@Override
	public double apply(Variable v) {
		return Math.sin(v.get());
	}
	
	@Override
	public MathFun _d(String varName) {
		return new FCos(getVarNames());
	}
}
