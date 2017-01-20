package edu.uta.futureye.core.intf;

import java.util.Map;

import edu.uta.futureye.function.intf.MathFunc;

public interface FiniteElement {
	MathFunc[] getShapeFunctions();
	int getNumberOfDOFs();
	Map<String, MathFunc> getCoordTransMap();
	String[] getArgsOrder();
	MathFunc getJacobian();
}
