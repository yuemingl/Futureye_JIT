package edu.uta.futureye.core.intf;

import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.function.intf.VectorMathFunc;

public interface LHSVecExpr {
	MathFunc apply(VectorMathFunc u, VectorMathFunc v);
}
