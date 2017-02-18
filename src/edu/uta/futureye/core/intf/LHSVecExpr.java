/**
 * Copyright (c) 2010, nkliuyueming@gmail.com. All rights reserved.
 * 
 * 
 */
package edu.uta.futureye.core.intf;

import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.function.intf.VecMathFunc;

/**
 * The left hand side (LHS) expression of a weak form with
 * vector valued shape functions
 *
 */
public interface LHSVecExpr {
	
	/**
	 * Return the bilinear form of the left hand side
	 * @param u
	 * @param v
	 * @return
	 */
	MathFunc apply(VecMathFunc u, VecMathFunc v);
}
