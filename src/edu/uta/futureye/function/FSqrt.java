package edu.uta.futureye.function;

import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.basic.FX;
import edu.uta.futureye.function.intf.MathFunc;

public class FSqrt extends AbstractSimpleMathFunc {

	public FSqrt(String funcName, String varName) {
		super("sqrt", varName);
	}

	@Override
	public double apply(Variable v) {
		return Math.sqrt(v.get(varName));
	}
	
	@Override
	public MathFunc diff(String varName) {
		return FC.c(0.5).M(pow(new FX(varName),-0.5)).M(f.diff(varName));
	}
	@Override
	public int getOpOrder() {
		return OP_ORDER1;
	}
	@Override
	public String toString() {
		return "sqrt("+f.toString()+")";
	}
	@Override
	public MathFunc copy() {
		// TODO Auto-generated method stub
		return null;
	}
	@Override
	public double apply(double... args) {
		// TODO Auto-generated method stub
		return args[argIdx];
	}

}
