package edu.uta.futureye.function.basic;

import edu.uta.futureye.function.AbstractFunction;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.util.Constant;

/**
 * f(x)=x
 * 
 * @author liuyueming
 *
 */
public class FX extends AbstractFunction{
	/**
	 * Use this to get f(x)=x, instead of construct a new FX, 
	 * it will be faster and memory saving :)
	 */
	public final static FX fx = new FX(Constant.x); 
	
	/**
	 * Use this to construct f(varName)
	 */
	public FX(String varName) {
		super(varName);
	}

	@Override
	public Function d(String varName) {
		if( this.varNames.contains(varName))
			return FC.c1;
		else
			return FC.c0;
	}

	@Override
	public double value(Variable v) {
		return v.get(varNames.get(0));
	}
	
	public String toString() {
		return varNames.get(0);
	}
}
