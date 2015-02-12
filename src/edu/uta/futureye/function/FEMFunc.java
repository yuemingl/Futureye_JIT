package edu.uta.futureye.function;

import symjava.bytecode.BytecodeFunc;
import symjava.numeric.NumInt;
import symjava.symbolic.Expr;
import symjava.symbolic.Func;
import symjava.symbolic.Int;
import symjava.symbolic.utils.Utils;

public class FEMFunc extends Func implements BytecodeFunc {
	NumInt integrate;
	
//	public FEMFunc(String name, Expr expr) {
//		super(name, expr);
//	}
//	
//	public FEMFunc(String name, Expr expr, Expr[] args) {
//		super(name, expr, args);
//	}
	
	public FEMFunc(Int integrate) {
		super(integrate.getLabel(), integrate);
		this.integrate = new NumInt(integrate);
	}

	@Override
	public double apply(double... args) {
		return integrate.eval(args);
	}

}
