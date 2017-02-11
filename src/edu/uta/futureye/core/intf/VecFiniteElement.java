package edu.uta.futureye.core.intf;

import java.util.Map;

import edu.uta.futureye.core.Element;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.function.intf.VectorMathFunc;

public interface VecFiniteElement {
	/**
	 * Return the number of Degree of Freedom (DOF) of this finite element
	 * @return
	 */
	int getNumberOfDOFs();
	
	/**
	 * Return all the expression of shape functions
	 * @return
	 */
	VectorMathFunc[] getShapeFunctions();
	
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
	//remove this function to assembler since the local to global index of DOF is better handled in assembler
	//void assignTo(Element e);
	
	/**
	 * Return the boundary FiniteElement object associated with the current VecFiniteElement object
	 * @return
	 */
	VecFiniteElement getBoundaryFE();
	
	/**
	 * check if the dot product of two vector valued shape functions equals 0 by
	 * providing the inedx of DOF.
	 * @param idx1
	 * @param idx2
	 * @return
	 */
	boolean isDOFCoupled(int idx1, int idx2);
}
