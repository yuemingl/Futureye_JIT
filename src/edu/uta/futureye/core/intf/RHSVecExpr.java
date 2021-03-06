/**
 * Copyright (c) 2010, nkliuyueming@gmail.com. All rights reserved.
 * 
 * 
 */
package edu.uta.futureye.core.intf;

import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.function.intf.VecMathFunc;

/**
 * The right hand side (RHS) expression of a weak form with
 * vector valued shape functions
 *
 */
public interface RHSVecExpr {
	
	/**
	 * Return the linear form of the right hand side
	 * @param v
	 * @return
	 */
	MathFunc apply(VecMathFunc v);
}
