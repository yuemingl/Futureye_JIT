package edu.uta.futureye.lib.element;


import java.util.Map;

import edu.uta.futureye.core.DOF;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.RectAreaCoord;
import edu.uta.futureye.core.Vertex;
import edu.uta.futureye.core.intf.FiniteElement;
import edu.uta.futureye.function.basic.FX;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.util.container.VertexList;

/**
 * Bilinear shape functions on a rectangle
 * 
 *  s
 *  ^
 *  |
 *  |
 * 
 *  4-----3
 *  |     |
 *  |     |
 *  1-----2  --> r
 * -1  0  1
 *
 * N1 = (1-r)*(1-s)/4
 * N2 = (1+r)*(1-s)/4
 * N3 = (1+r)*(1+s)/4
 * N4 = (1-r)*(1+s)/4
 * 
*/
public class FEBilinearRectangle implements FiniteElement {
	RectAreaCoord coord;
	
	//Construct a function with the coordinate of points in an element as parameters
	String[] argsOrder;
	
	public int nDOFs = 4;
	MathFunc[] shapeFuncs = new MathFunc[nDOFs];

	public FEBilinearRectangle() {
		FX x1 = new FX("x1");
		FX x2 = new FX("x2");
		FX x3 = new FX("x3");
		FX x4 = new FX("x4");
		FX y1 = new FX("y1");
		FX y2 = new FX("y2");
		FX y3 = new FX("y3");
		FX y4 = new FX("y4");

		this.coord = new RectAreaCoord(x1, x2, x3, x4, y1, y2, y3, y4);
		MathFunc r = coord.getCoordR();
		MathFunc s = coord.getCoordS();
		
		this.argsOrder = new String[]{x1,x2,x3,x4,y1,y2,y3,y4,r,s};

		this.shapeFuncs[0] = (1-r)*(1-s)/4;
		this.shapeFuncs[1] = (1+r)*(1-s)/4;
		this.shapeFuncs[2] = (1+r)*(1+s)/4;
		this.shapeFuncs[3] = (1-r)*(1+s)/4;
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
	public Map<String, MathFunc> getCoordTransMap() {
		return this.coord.getCoordTransMap();
	}

	@Override
	public String[] getArgsOrder() {
		return this.argsOrder;
	}
	
	@Override
	public MathFunc getJacobian() {
		return this.coord.getJacobian();
	}

	public void assignTo(Element e) {
		e.clearAllDOF();
		VertexList vertices = e.vertices();
		for(int j=1;j<=vertices.size();j++) {
			Vertex v = vertices.at(j);
			//Assign shape function to DOF
			DOF dof = new DOF(
						j, //Local DOF index
						v.globalNode().getIndex(), //Global DOF index, take global node index
						null //Shape function is no longer used?  
						);
			e.addNodeDOF(j, dof);
		}
	}

	@Override
	public FiniteElement getBoundaryFE() {
		return new FELinearLine2D();
	}

}
