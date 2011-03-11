package edu.uta.futureye.function.basic;

import edu.uta.futureye.function.AbstractFunction;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.intf.Function;

public class FLinear1D  extends AbstractFunction{
	protected double x1,x2,y1,y2;
	
	public FLinear1D(double x1, double y1, 
			double x2, double y2) {
		this.x1 = x1;
		this.y1 = y1;

		this.x2 = x2;
		this.y2 = y2;
	}
	
	@Override
	public double value(Variable v) {
		double x = v.get(varNames().get(0));
		return (x-x1)*(y2-y1)/(x2-x1) + y1;
	}

	@Override
	public Function d(String varName) {
		if(this.varNames().contains(varName))
			return new FC((y2-y1)/(x2-x1));
		else
			return new FC(0.0);
	}
	
	public String toString() {
		return (y2-y1)/(x2-x1)+"*("+varNames().get(0)+"-"+x1+")+"+y1;
	}
}
