package edu.uta.futureye.core.intf;

import edu.uta.futureye.function.intf.MathFunc;

public interface LHSExpr {
	MathFunc apply(MathFunc u, MathFunc v);
}
