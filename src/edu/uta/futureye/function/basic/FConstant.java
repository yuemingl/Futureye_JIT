package edu.uta.futureye.function.basic;

import java.util.LinkedList;
import java.util.List;

import edu.uta.futureye.function.DerivativeIndicator;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.intf.FunctionDerivable;

public class FConstant extends FAbstract{
	protected double val = 0.0;
	
	public FConstant(double val) {
		List<String> varNames = new LinkedList<String>();
		//TODO !!!
		setVarNames(varNames);
		this.val = val;
	}
	
	@Override
	public double value(Variable v) {
		return val;
	}

	public double value() {
		return val;
	}
	
	@Override
	public FunctionDerivable derivative(DerivativeIndicator di) {
		return new FConstant(0.0);
	}
	
	public String toString() {
		return String.valueOf(val);
	}
}
