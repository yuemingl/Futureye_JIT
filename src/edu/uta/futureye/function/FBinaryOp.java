package edu.uta.futureye.function;

import edu.uta.futureye.function.intf.MathFunc;

public abstract class FBinaryOp extends AbstractMathFunc {
	public MathFunc arg1;
	public MathFunc arg2;
	
	public FBinaryOp(MathFunc arg1, MathFunc arg2) {
		this.arg1 = arg1;
		this.arg2 = arg2;
	}
	
	public MathFunc[] args() {
		return new MathFunc[] { arg1, arg2 };
	}
	
	public MathFunc lhs() {
		return arg1;
	}
	
	public MathFunc rhs() {
		return arg2;
	}
}
