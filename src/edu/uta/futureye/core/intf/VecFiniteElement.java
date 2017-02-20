/**
 * Copyright (c) 2010, nkliuyueming@gmail.com. All rights reserved.
 * 
 * 
 */
package edu.uta.futureye.core.intf;

import edu.uta.futureye.core.DOF;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.NodeType;
import edu.uta.futureye.core.geometry.GeoEntity;
import edu.uta.futureye.function.intf.VecMathFunc;

/**
 * Finite element with vector valued shape functions
 *
 */
public interface VecFiniteElement {
	/**
	 * Return the number of Degree of Freedoms (DOFs) of this finite element
	 * @return
	 */
	int getNumberOfDOFs();
	
	/**
	 * Return all the expression of shape functions
	 * @return
	 */
	VecMathFunc[] getShapeFunctions();
	
	/**
	 * Return the geometry entity on element e giving the local index of shape function
	 * @param localIndex
	 * @return
	 */
	NodeType getDOFType(Element e, int localIndex);

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
	 * Return the index of component of the vector valued function by giving 
	 * a local index of a shape function
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
	int getGlobalIndex(Mesh mesh, Element e, int localIndex);
	
	/**
	 * Get all the degree of freedoms (DOF) on a given mesh for this finite element
	 * @param mesh
	 * @return
	 */
	int getTotalNumberOfDOFs(Mesh mesh);
	
	/**
	 * Get the degree of freedoms (DOF) for the given vector valued function (VVF) component
	 * on a given mesh for this finite element. For example,
	 * 1 - the degree of freedoms for the velocity component v
	 * 2 - the degree of freedoms for the velocity component u
	 * 3 - the degree of freedoms for the pressure component p
	 * @param mesh
	 * @param nVVFComponentIndex
	 * @return
	 */
	int getNumberOfNOFs(Mesh mesh, int nVVFComponentIndex);
	
	/**
	 * Return the coordinate transformation object used in this finite element
	 * by giving the index of component of the vector finite element
	 * @return
	 */
	//CoordTrans getCoordTrans(int index);
	CoordTrans getCoordTrans();
	
	/**
	 * Return the ODF object at the given local index
	 * @param localIndex
	 * @return
	 */
	DOF getDOF(int localIndex);
	
	/**
	 * Return the geometry entity at the given local index
	 * @param e
	 * @param localIndex
	 * @return
	 */
	GeoEntity getGeoEntity(Element e, int localIndex);
}
