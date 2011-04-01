package edu.uta.futureye.lib.element;

import edu.uta.futureye.core.Element;

public interface FiniteElementType {
	void assign(Element e);
	int getVectorShapeFunctionDim();
	int getDOFNumOnElement(int vsfDim);
}
