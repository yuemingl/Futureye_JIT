package edu.uta.futureye.function.basic;

import edu.uta.futureye.function.AbstractFunction;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.util.Constant;

public class FAx extends AbstractFunction {
	protected double a;
	
	public FAx(double a) {
		varNames.add(Constant.x);
		this.a = a;
	}
	
	public FAx(String varName, double a) {
		super(varName);
		this.a = a;
	}
	
	@Override
	public Function _d(String varName) {
		if(this.varNames().contains(varName))
			return new FC(a);
		else
			return FC.c0;
	}

	@Override
	public double value(Variable v) {
		return a*v.get(varNames().get(0));
	}
	
	@Override
	public Function copy() {
		return new FAx(this.varNames.get(0),a);
	}
	
	@Override
	public String toString() {
		if(Double.compare(a, 1.0) == 0)
			return " "+varNames().get(0)+" ";
		else if(Double.compare(a, 0.0) == 0)
			return " 0.0 ";
		return " "+a+"*"+varNames().get(0)+" ";
	}
	
}
