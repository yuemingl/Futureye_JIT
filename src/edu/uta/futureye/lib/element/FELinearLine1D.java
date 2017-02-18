/**
 * Copyright (c) 2010, nkliuyueming@gmail.com. All rights reserved.
 * 
 * 
 */
package edu.uta.futureye.lib.element;

import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Line1DCoord;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.intf.CoordTrans;
import edu.uta.futureye.core.intf.FiniteElement;
import edu.uta.futureye.function.basic.FX;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.util.container.VertexList;

/**
 * Linear line element in a 1D space
 * 
 * shape functions
 * 
 *  1-----2  -->r
 * -1  0  1
 * 
 * N1 = (1-r)/2
 * N2 = (1+r)/2
 * 
 */
public class FELinearLine1D implements FiniteElement {
	Line1DCoord coord;
	
	//Construct a function with the coordinate of points in an element as parameters
	String[] argsOrder;
	
	public int nDOFs = 2;
	MathFunc[] shapeFuncs = new MathFunc[nDOFs];

	public FELinearLine1D() {
		FX x1 = new FX("x1");
		FX x2 = new FX("x2");
		
		this.coord = new Line1DCoord(x1, x2);
		MathFunc r = coord.getCoordR();
		
		argsOrder = new String[]{x1,x2,r};
		
		shapeFuncs[0] = 0.5*(1-r);
		shapeFuncs[1] = 0.5*(1+r);
	}

	@Override
	public MathFunc[] getShapeFunctions() {
		return this.shapeFuncs;
	}

	@Override
	public int getNumberOfDOFs() {
		return this.nDOFs;
	}

	@Override
	public String[] getArgsOrder() {
		return this.argsOrder;
	}

	@Override
	public FiniteElement getBoundaryFE() {
		throw new UnsupportedOperationException();
	}

	@Override
	public int getGlobalIndex(Mesh mesh, Element e, int localIndex) {
		VertexList vertices = e.vertices();
		return vertices.at(localIndex).globalNode().globalIndex;
	}

	@Override
	public int getTotalNumberOfDOFs(Mesh mesh) {
		return mesh.getNodeList().size();
	}

	@Override
	public CoordTrans getCoordTrans() {
		return this.coord;
	}
}
