package edu.uta.futureye.bytecode;

import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.function.intf.MathFunc;

/**
 * Represents a compiled function from an object of MathFunc
 *
 */
public abstract class CompiledFunc {
	protected MathFunc[] funcRefs;

	/**
	 * This function should be implemented during compilation
	 * @param e
	 * @param n
	 * @param args
	 * @return
	 */
	public abstract double apply(Element e, Node n, double ...args);
	
	/**
	 * A simplified method for applying a function without 
	 * parameters Element and Node
	 * @param args
	 * @return
	 */
	public double apply(double ...args) {
		return apply(null, null, args);
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