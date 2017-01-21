package edu.uta.futureye.core.intf;

import java.util.Map;

import edu.uta.futureye.core.Element;
import edu.uta.futureye.function.intf.MathFunc;

public interface FiniteElement {
	/**
	 * Return the number of Degree of Freedom (DOF) of this finite element
	 * @return
	 */
	int getNumberOfDOFs();
	
	/**
	 * Return all the expression of shape functions
	 * @return
	 */
	MathFunc[] getShapeFunctions();
	
	/**
	 * Return the coordinate transformation map between physical coordinate and reference coordinate
	 * @return
	 */
	Map<String, MathFunc> getCoordTransMap();
	
	/**
	 * Return the order of all the arguments
	 * @return
	 */
	String[] getArgsOrder();
	
	/**
	 * Return the Jacobian expression
	 * @return
	 */
	MathFunc getJacobian();
	
	/**
	 * Associate this FiniteElement object to an Element object which contains geometry information
	 * 
	 * @param e
	 */
	void assignTo(Element e);
	
	/**
	 * Return the boundary FiniteElement object associated with the current FiniteElement object
	 * @return
	 */
	FiniteElement getBoundaryFE();
}
