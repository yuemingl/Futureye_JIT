package edu.uta.futureye.core.intf;

import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.function.intf.VecMathFunc;

public interface RHSVecExpr {
	MathFunc apply(VecMathFunc v);
}
