package edu.uta.futureye.function.intf;

import edu.uta.futureye.core.Element;

public interface ShapeFunction extends FunctionDerivable {
	public void asignElement(Element e);
	public ShapeFunction restrictTo(int funIndex);
}
