package edu.uta.futureye.function.basic;

import edu.uta.futureye.function.DerivativeIndicator;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.intf.FunctionDerivable;

public class FAx extends FAbstract {
	protected double a = 1.0;

	public FAx() {
	}
	
	public FAx(String varName, double coef) {
		super(varName);
		this.a = coef;
	}
	
	public FAx(double coef) {
		this.a = coef;
	}
	
	@Override
	public FunctionDerivable derivative(DerivativeIndicator di) {
		Integer degree = di.get().get(varNames().get(0));
		if( degree != null && degree == 1)
			return new FConstant(a);
		else
			return new FConstant(0.0);
	}

	@Override
	public double value(Variable v) {
		return a*v.get(varNames().get(0));
	}
	
	public String toString() {
		if(Double.compare(a, 1.0) == 0)
			return " "+varNames().get(0)+" ";
		else if(Double.compare(a, 0.0) == 0)
			return " 0.0 ";
		return " "+a+"*"+varNames().get(0)+" ";
	}
	
}
