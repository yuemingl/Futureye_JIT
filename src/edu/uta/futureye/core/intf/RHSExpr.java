package edu.uta.futureye.core.intf;

import edu.uta.futureye.function.intf.MathFunc;

public interface RHSExpr {
	MathFunc apply(MathFunc v);
}
