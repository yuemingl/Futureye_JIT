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
	 * Used to form f(x)=x, instead of construct a new FX object, 
	 * it will be faster and memory saving :)
	 */
	public final static FX fx = new FX(Constant.x); 
	
	/**
	 * Different variable names
	 */
	public final static FX fy = new FX(Constant.y); 
	public final static FX fz = new FX(Constant.z); 
	
	public final static FX fr = new FX(Constant.r); 
	public final static FX fs = new FX(Constant.s); 
	public final static FX ft = new FX(Constant.t); 
	
	/**
	 * Use this to construct f(varName)
	 */
	public FX(String varName) {
		super(varName);
	}

	@Override
	public Function _d(String varName) {
		if( this.varNames.contains(varName))
			return FC.c1;
		else
			return FC.c0;
	}

	@Override
	public double value(Variable v) {
		return v.get(varNames.get(0));
	}
	
	@Override
	public int getOpOrder() {
		return OP_ORDER0;
	}
	
	@Override
	public Function copy() {
		return new FX(this.varNames.get(0));
	}
	
	@Override
	public String toString() {
		return varNames.get(0);
	}
}
