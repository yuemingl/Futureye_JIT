package edu.uta.futureye.function.basic;

import edu.uta.futureye.function.DerivativeIndicator;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.intf.FunctionDerivable;

public class FX extends FAbstract{

	public FX() {
	}
	
	public FX(String varName) {
		super(varName);
	}

	
	@Override
	public FunctionDerivable derivative(DerivativeIndicator di) {
		Integer degree = di.get().get(varNames().get(0));
		if( degree != null && degree == 1)
			return new FConstant(1.0);
		else
			return new FConstant(0.0);
	}

	@Override
	public double value(Variable v) {
		return v.get(varNames().get(0));
	}
	
	public String toString() {
		return varNames().get(0);
	}
}
