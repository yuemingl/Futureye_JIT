package edu.uta.futureye.function;

import edu.uta.futureye.function.intf.MathFun;

public abstract class FBinaryOp extends AbstractMathFun {
	public MathFun arg1;
	public MathFun arg2;
	
	public FBinaryOp(MathFun arg1, MathFun arg2) {
		this.arg1 = arg1;
		this.arg2 = arg2;
	}
	
	public MathFun[] args() {
		return new MathFun[] { arg1, arg2 };
	}
	
	public MathFun lhs() {
		return arg1;
	}
	
	public MathFun rhs() {
		return arg2;
	}
}
