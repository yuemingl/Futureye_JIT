package edu.uta.futureye.function.basic;

import edu.uta.futureye.function.AbstractFunction;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.intf.Function;

public class FAxpb extends AbstractFunction {
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
	public Function d(String varName) {
		if(this.varNames().contains(varName))
			return new FC(a);
		else
			return FC.c0;
	}

	@Override
	public double value(Variable v) {
		double rlt = a*v.get(varNames.get(0))+b;
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
