package edu.uta.futureye.function.basic;

import java.util.List;

import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.function.AbstractMathFunc;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.intf.MathFunc;

public class FCos extends AbstractMathFunc{
	
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
		return Math.cos(v.get());
	}

	@Override
	public double apply(Element e, Node n, double... args) {
		return Math.cos(args[0]);
	}

	@Override
	public double apply(double... args) {
		return apply(null, null, args);
	}
	
	@Override
	public MathFunc diff(String varName) {
		return (new FSin(getVarNames())).M(-1);
	}

}
