/**
 * Copyright (c) 2010, nkliuyueming@gmail.com. All rights reserved.
 * 
 * 
 */
package edu.uta.futureye.core.intf;

import edu.uta.futureye.function.intf.MathFunc;

/**
 * The right hand side (RHS) expression of a weak form
 *
 */
public interface RHSExpr {

	/**
	 * Return the linear form of the right hand side of a weak form
	 * 
	 * @param v
	 * @return
	 */
	MathFunc apply(MathFunc v);
}
