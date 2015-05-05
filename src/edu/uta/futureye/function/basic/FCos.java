package edu.uta.futureye.function.basic;

import java.util.List;

import edu.uta.futureye.function.AbstractMathFun;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.intf.MathFunc;

public class FCos extends AbstractMathFun{
	
	public FCos() {
		super("x");
	}
	
	public FCos(String varName) {
		super(varName);
	}
	
	public FCos(List<String> varNames) {
		super(varNames);
	}
	
	@Override
	public double apply(Variable v) {
		return Math.sin(v.get());
	}
	
	@Override
	public MathFunc _d(String varName) {
		return (new FSin(getVarNames())).M(-1);
	}

}
