package edu.uta.futureye.function.basic;

import edu.uta.futureye.function.AbstractSimpleMathFunc;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.util.Constant;

/**
 * f(x) = sin(x)
 *
 */
public class FSin extends AbstractSimpleMathFunc {
	public FSin() {
		super("sin", Constant.x);
	}
	
	public FSin(String varName) {
		super("sin", varName);
	}
	
	@Override
	public double apply(Variable v) {
		return Math.sin(v.get());
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
	public String getExpr() {
		return "sin("+varName+")";
	}
	
	@Override
	public String toString() {
		return "sin("+varName+")";
	}

}
