/**
 * Copyright (c) 2010, nkliuyueming@gmail.com. All rights reserved.
 * 
 * 
 */
package edu.uta.futureye.bytecode;

import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.lib.assembler.AssembleParam;

/**
 * This class represents a compiled function from an object of MathFunc
 *
 */
public abstract class CompiledFunc {
	//references to MathFunc objects
	protected MathFunc[] funcRefs;

	/**
	 * This function should be implemented during compilation
	 * @param ap - This parameter is supposed to be passed into the function during assembling process
	 * @param args - arguments defined in MathFunc
	 * @return
	 */
	public abstract double apply(AssembleParam ap, double ...args);
	
	public double call(AssembleParam ap, double ...args) {
		return apply(ap, args);
	}

	/**
	 * A simplified method for applying a function without 
	 * parameters AssembleParam
	 * @param args
	 * @return
	 */
	public double apply(double ...args) {
		return apply(null, args);
	}
	
	public double call(double ...args) {
		return apply(null, args);
	}

	/**
	 * Set the references to functions before compilation for
	 * further calling from apply(...)
	 * @param funcs
	 */
	public void setFuncRefs(MathFunc[] funcs) {
		this.funcRefs = funcs;
	}
}