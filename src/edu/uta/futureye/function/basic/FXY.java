package edu.uta.futureye.function.basic;

import java.util.List;

import edu.uta.futureye.function.AbstractFunction;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.util.Constant;

/**
 * f(x,y) = c1*x + c2*y + c3
 * 
 * @author liuyueming
 *
 */
public class FXY extends AbstractFunction{
	protected double c1,c2,c3=0.0;

	public FXY(double c1,double c2) {
		varNames.add(Constant.x);
		varNames.add(Constant.y);
		this.c1 = c1;
		this.c2 = c2;
	}
	
	public FXY(double c1,double c2,double c3) {
		varNames.add(Constant.x);
		varNames.add(Constant.y);
		this.c1 = c1;
		this.c2 = c2;
		this.c3 = c3;
	}	
	
	public FXY(List<String> varNames,double c1,double c2) {
		super(varNames);
		this.c1 = c1;
		this.c2 = c2;
	}
	
	public FXY(List<String> varNames,double c1,double c2,double c3) {
		super(varNames);
		this.c1 = c1;
		this.c2 = c2;
		this.c3 = c3;
	}	
	
	@Override
	public Function _d(String varName) {
		if(varNames.get(0).equals(varName))
			return new FC(c1);
		else if(varNames.get(1).equals(varName)) {
			return new FC(c2);
		}
		return FC.c0;
	}

	@Override
	public double value(Variable v) {
		return c1 * v.get(varNames.get(0)) + c2 * v.get(varNames.get(1)) + c3;
	}
	
	public String toString() {
		String s1 = "";
		if(Math.abs(c1) > Constant.eps)
			s1 = s1 + "*";
		String s2 = "";
		if(Math.abs(c2) > Constant.eps)
			s2 = s2 + "*";
		if(s1.isEmpty() && !s2.isEmpty())
			return s2+varNames.get(1);
		else if(!s1.isEmpty() && s2.isEmpty())
			return s1+varNames.get(0);
		else
			return s1+varNames.get(0)+" + " +s2+varNames.get(1);
		//TODO c3
	}
}