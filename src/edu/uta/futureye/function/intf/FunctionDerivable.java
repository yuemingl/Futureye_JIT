package edu.uta.futureye.function.intf;

import edu.uta.futureye.function.DerivativeIndicator;

public interface FunctionDerivable extends Function {
	public FunctionDerivable derivative(DerivativeIndicator di);
}
