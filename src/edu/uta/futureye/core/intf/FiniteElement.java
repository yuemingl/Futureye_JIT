/**
 * Copyright (c) 2010, nkliuyueming@gmail.com. All rights reserved.
 * 
 * 
 */
package edu.uta.futureye.core.intf;

import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.function.intf.MathFunc;

/**
 * Finite element interface which provides the information about
 * the degree of freedoms, shape functions, local-to-global index and 
 * boundary elements.
 * 
 * A user defined finite element should implements the declared
 * functions in this interface.
 *
 */
public interface FiniteElement {
	/**
	 * Return the number of degree of freedoms (DOFs) of this finite element
	 * TODO: change to getDegreeOfFreedom() ???
	 * @return
	 */
	int getNumberOfDOFs();

	/**
	 * Return all the shape functions of this finite element
	 * @return
	 */
	MathFunc[] getShapeFunctions();

	/**
	 * Return the order of all the arguments in a string array
	 * @return
	 */
	String[] getArgsOrder();

	/**
	 * Get the global index of a degree of freedom (DOF) by giving the local index
	 * @param mesh
	 * @param e
	 * @param localIndex
	 * @return
	 */
	int getGlobalIndex(Mesh mesh, Element e, int localIndex);

	/**
	 * Get the total number of degree of freedoms (DOFs) on a given mesh for this finite element
	 * @param mesh
	 * @return
	 */
	int getTotalNumberOfDOFs(Mesh mesh);

	/**
	 * Returns the boundary FiniteElement object associated with this FiniteElement object.
	 * For example, a boundary finite element is an element on a line for a 2D element.
	 * For a 3D element, the boundary finite element is an element on a 2D face.
	 * @return
	 */
	FiniteElement getBoundaryFE();
	
	/**
	 * Return the coordinate transformation object used in this finite element
	 * @return
	 */
	CoordTrans getCoordTrans();
}
