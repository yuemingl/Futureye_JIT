package edu.uta.futureye.function.basic;

import edu.uta.futureye.function.AbstractSimpleMathFunc;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.util.Constant;

/**
 * f(x) = tan(x)
 *
 */
public class FTan extends AbstractSimpleMathFunc {

	public FTan() {
		super("tan", Constant.x);
	}
	
	public FTan(String varName) {
		super("tan", varName);
	}

//	@Override
//	public MathFunc copy() {
//		FTan ret = new FTan(this.varName);
//		ret.argIdx = this.argIdx;
//		return ret;
//	}

	@Override
	public double apply(double... args) {
		return Math.tan(args[argIdx]);
	}

	@Override
	public double apply(Variable v) {
		return Math.tan(v.get(varName));
	}
	
	@Override
	public MathFunc diff(String varName) {
		//1 + tan^2(x) 
		if(varName.equals(this.varName))
			return this.M(this).A(1);
		else
			return FC.C0;
	}
}
