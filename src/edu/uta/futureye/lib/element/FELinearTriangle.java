package edu.uta.futureye.lib.element;


import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.TriAreaCoord;
import edu.uta.futureye.core.intf.CoordTrans;
import edu.uta.futureye.core.intf.FiniteElement;
import edu.uta.futureye.function.basic.FX;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.util.container.VertexList;

/**
 * Linear finite element on a triangle element.
 *
 */
public class FELinearTriangle implements FiniteElement {
	TriAreaCoord coord;
	
	//Construct a function with the coordinate of points in an element as parameters
	String[] argsOrder;
	
	public int nDOFs = 3;
	MathFunc[] shapeFuncs = new MathFunc[nDOFs];

	public FELinearTriangle() {
		FX x1 = new FX("x1");
		FX x2 = new FX("x2");
		FX x3 = new FX("x3");
		FX y1 = new FX("y1");
		FX y2 = new FX("y2");
		FX y3 = new FX("y3");
		
		this.coord = new TriAreaCoord(x1, x2, x3, y1, y2, y3);
		MathFunc r = coord.getCoordR();
		MathFunc s = coord.getCoordS();
		
		this.argsOrder = new String[]{x1, x2, x3, y1, y2, y3, r, s, "t"};
		
		//shape functions
		this.shapeFuncs[0] = r;
		this.shapeFuncs[1] = s;
		this.shapeFuncs[2] = 1 - r - s;
	}

	@Override
	public int getNumberOfDOFs() {
		return this.nDOFs;
	}

	@Override
	public MathFunc[] getShapeFunctions() {
		return this.shapeFuncs;
	}

	@Override
	public String[] getArgsOrder() {
		return this.argsOrder;
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
	public FiniteElement getBoundaryFE() {
		return new FELinearLine2D();
	}

	@Override
	public CoordTrans getCoordTrans() {
		return this.coord;
	}
}
