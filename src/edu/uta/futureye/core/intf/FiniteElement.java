package edu.uta.futureye.core.intf;

import java.util.Map;

import edu.uta.futureye.core.Element;
import edu.uta.futureye.function.intf.MathFunc;

public interface FiniteElement {
	int getNumberOfDOFs();
	
	MathFunc[] getShapeFunctions();
	
	Map<String, MathFunc> getCoordTransMap();
	
	String[] getArgsOrder();
	
	MathFunc getJacobian();
	
	void assignTo(Element e);
}
