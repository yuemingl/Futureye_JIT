/**
 * Copyright (c) 2010, nkliuyueming@gmail.com. All rights reserved.
 * 
 * 
 */
package edu.uta.futureye.core.intf;

import edu.uta.futureye.function.intf.MathFunc;

/**
 * The left hand side (LHS) expression of a weak form
 *
 */
public interface LHSExpr {
	
	/**
	 * Return the bilinear form of the left hand side of a weak form
	 * 
	 * @param u
	 * @param v
	 * @return
	 */
	MathFunc apply(MathFunc u, MathFunc v);
}
