package edu.uta.futureye.function.basic;

import java.util.List;

import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.function.AbstractMathFun;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.intf.MathFunc;

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
	public double apply(Element e, Node n, double... args) {
		return Math.sin(args[0]);
	}

	@Override
	public double apply(double... args) {
		return apply(null, null, args);
	}
	
	@Override
	public MathFunc _d(String varName) {
		return new FCos(getVarNames());
	}
}
