package edu.uta.futureye.function.basic;

import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.function.AbstractSimpleMathFunc;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.intf.MathFunc;

public class FSin extends AbstractSimpleMathFunc {
	
	public FSin() {
		super("sin", "x");
	}
	
	public FSin(String varName) {
		super("sin", varName);
	}
	
	@Override
	public double apply(Variable v) {
		return Math.sin(v.get());
	}

	@Override
	public double apply(Element e, Node n, double... args) {
		return Math.sin(args[argIdx]);
	}

	@Override
	public double apply(double... args) {
		return Math.sin(args[argIdx]);
	}
	
	@Override
	public MathFunc diff(String varName) {
		if(varName.equals(this.varName))
			return new FCos(this.varName);
		else
			return FC.C0;
	}

	@Override
	public MathFunc copy() {
		FSin ret = new FSin(this.varName);
		ret.argIdx = this.argIdx;
		return ret;
	}
}
