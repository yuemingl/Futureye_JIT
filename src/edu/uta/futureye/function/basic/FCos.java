package edu.uta.futureye.function.basic;

import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.function.AbstractSimpleMathFunc;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.intf.MathFunc;

public class FCos extends AbstractSimpleMathFunc{
	
	public FCos() {
		super("cos", "x");
	}
	
	public FCos(String varName) {
		super("cos", varName);
	}
	
	@Override
	public double apply(Variable v) {
		return Math.cos(v.get());
	}

	@Override
	public double apply(Element e, Node n, double... args) {
		return Math.cos(args[argIdx]);
	}

	@Override
	public double apply(double... args) {
		return apply(null, null, args);
	}
	
	@Override
	public MathFunc diff(String varName) {
		if(varName.equals(this.varName))
			return (new FSin(this.varName)).M(-1);
		else
			return FC.C0;
	}

	@Override
	public MathFunc copy() {
		FCos ret = new FCos(this.varName);
		ret.argIdx = this.argIdx;
		return ret;
	}

}
