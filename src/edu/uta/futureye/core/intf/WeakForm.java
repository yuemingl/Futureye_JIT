package edu.uta.futureye.core.intf;


import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.function.intf.ShapeFunction;

public interface WeakForm {
	static enum ItemType {Domain, Border};
	
	//////////////////////Common Approach/////////////////////////
	void setShapeFunction(ShapeFunction trial, int trialDofLocalIndex,
			ShapeFunction test, int testDofLocalIndex);
	
	Function leftHandSide(Element e, ItemType itemType);
	Function rightHandSide(Element e, ItemType itemType);
	
	//////////////////////Fast Approach//////////////////////////
	/**
	 * Assemble element e here, instead of providing left hand side
	 * and right hand side.
	 * 
	 * @param e
	 * @param globalStiff (I/O): Global stiff matrix 
	 * @param globalLoad (I/O): Global load vector
	 *   
	 */
	void assembleElement(Element e, 
			Matrix globalStiff, Vector globalLoad);
	
	///////////////////////////////////////////////////////////////
	
	/**
	 * Integrate on element
	 * 
	 * @param e
	 * @param fun: LHS or RHS
	 * @return
	 */
	double integrate(Element e, Function fun);
	
}
