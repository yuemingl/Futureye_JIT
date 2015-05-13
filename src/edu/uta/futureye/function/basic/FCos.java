package edu.uta.futureye.function.basic;

import edu.uta.futureye.function.AbstractSimpleMathFunc;
import edu.uta.futureye.function.FMath;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.util.Constant;

/**
 * f(x) = cos(x)
 *
 */
public class FCos extends AbstractSimpleMathFunc{
	
	public FCos() {
		super("cos", Constant.x);
	}
	
	public FCos(String varName) {
		super("cos", varName);
	}

	@Override
	public double apply(double... args) {
		return Math.cos(args[argIdx]);
	}
	
	@Override
	public MathFunc diff(String varName) {
		if(varName.equals(this.varName))
			return (new FSin(this.varName)).M(-1);
		else
			return FMath.C0;
	}

	@Override
	public String getExpr() {
		return "cos("+varName+")";
	}
	
	@Override
	public String toString() {
		return "cos("+varName+")";
	}

}
