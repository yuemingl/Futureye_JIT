package edu.uta.futureye.lib.element;

import edu.uta.futureye.core.Element;

public interface FiniteElementType {
	/**
	 * Associate degrees of freedom (DOF) to element e
	 * @param e
	 */
	void assignTo(Element e);
	
	int getVectorShapeFunctionDim();
	
	int getDOFNumOnElement(int vsfDim);
}
