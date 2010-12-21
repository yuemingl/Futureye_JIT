package edu.uta.futureye.core.intf;

import java.util.List;

import edu.uta.futureye.algebra.Matrix;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.function.intf.ShapeFunction;
import edu.uta.futureye.util.PairElementMatrix;

public interface WeakForm {
	public static enum ItemType {Domain, Border};
	
	public void setShapeFunction(ShapeFunction trial, int trialDofLocalIndex,
			ShapeFunction test, int testDofLocalIndex);
	
	public Function leftHandSide(Element e, ItemType itemType);
	public Function rightHandSide(Element e, ItemType itemType);
	
	/**
	 * Associate an element to the weak form. 
	 * The local matrix can be assembled in this step if necessary.
	 * (including element e itself and the border elements of e etc.)
	 * @param e
	 * @return 
	 */
	public List<PairElementMatrix> associateElement(Element e);
	
}
