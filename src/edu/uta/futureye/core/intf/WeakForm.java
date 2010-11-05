package edu.uta.futureye.core.intf;

import edu.uta.futureye.core.Element;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.function.intf.ShapeFunction;

public interface WeakForm {
	public static enum ItemType {Domain, Border};
	public void setShapeFunction(ShapeFunction trial, int trialDofLocalIndex,
			ShapeFunction test, int testDofLocalIndex);
	public Function leftHandSide(Element e, ItemType itemType);
	public Function rightHandSide(Element e, ItemType itemType);
	
}
