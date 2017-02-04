package edu.uta.futureye.core.intf;

import edu.uta.futureye.core.Element;
import edu.uta.futureye.function.intf.MathFunc;

public interface RHSExprWithElement {
	MathFunc apply(MathFunc v, Element e);
}
