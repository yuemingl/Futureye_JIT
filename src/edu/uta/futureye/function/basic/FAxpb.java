package edu.uta.futureye.function.basic;

import edu.uta.futureye.function.DerivativeIndicator;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.intf.FunctionDerivable;

public class FAxpb extends FAbstract {
	protected double a;
	protected double b;

	public FAxpb(double a, double b) {
		this.a = a;
		this.b = b;
	}
	
	public FAxpb(String varName, double a, double b) {
		super(varName);
		this.a = a;
		this.b = b;
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
		double rlt = a*v.get(varNames().get(0))+b;
		return rlt;
	}
	
	public String toString() {
		if(Double.compare(a, 1.0) == 0) {
			if(Double.compare(b, 0.0) == 0)
				return " "+varNames().get(0)+" ";
			else
				return " "+varNames().get(0)+"+"+b+" ";
		} else if(Double.compare(a, 0.0) == 0) {
			if(Double.compare(b, 0.0) == 0)
				return " 0.0 ";
			else
				return " "+b+" ";			
		}
		return "("+a+"*"+varNames().get(0)+"+"+b+")";
	}
}
