package edu.uta.futureye.core.intf;

import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.function.intf.VecMathFunc;

/**
 * Finite element with vector valued shape functions
 *  
 */
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
	VecMathFunc[] getShapeFunctions();
		
	/**
	 * Return the order of all the arguments
	 * @return
	 */
	String[] getArgsOrder();
		
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
	
	/**
	 * Return the component index of the vector valued function by giving 
	 * a local index of shape functions
	 * @param localIndex
	 * @return
	 */
	int getVVFComponentIndex(int localIndex);
	
	/**
	 * Get the global index of a DEF
	 * @param mesh
	 * @param e
	 * @param localIndex
	 * @return
	 */
	public int getGlobalIndex(Mesh mesh, Element e, int localIndex);
	
	/**
	 * Get the total number of DOFs on a given mesh for this finite element
	 * @param mesh
	 * @return
	 */
	public int getTotalNumberOfDOFs(Mesh mesh);
	
	/**
	 * Return the coordinate transformation object used in this finite element
	 * by giving the index of component of the vector finite element
	 * @return
	 */
	//CoordTrans getCoordTrans(int index);
	CoordTrans getCoordTrans();

}
