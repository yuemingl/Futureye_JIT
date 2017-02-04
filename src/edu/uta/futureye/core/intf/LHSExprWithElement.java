package edu.uta.futureye.core.intf;

import edu.uta.futureye.core.Element;
import edu.uta.futureye.function.intf.MathFunc;

public interface LHSExprWithElement {
	MathFunc apply(MathFunc u, MathFunc v, Element e);

}
