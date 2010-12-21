package edu.uta.futureye.core;

import java.util.List;

import edu.uta.futureye.algebra.Matrix;
import edu.uta.futureye.core.intf.WeakForm;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.function.intf.ShapeFunction;
import edu.uta.futureye.util.PairElementMatrix;

public class AbstractWeakForm implements WeakForm {

	@Override
	public List<PairElementMatrix> associateElement(Element e) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Function leftHandSide(Element e, ItemType itemType) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Function rightHandSide(Element e, ItemType itemType) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void setShapeFunction(ShapeFunction trial, int trialDofLocalIndex,
			ShapeFunction test, int testDofLocalIndex) {
		// TODO Auto-generated method stub

	}

}
