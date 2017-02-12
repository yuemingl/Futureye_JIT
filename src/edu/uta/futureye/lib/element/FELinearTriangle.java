package edu.uta.futureye.lib.element;


import java.util.Map;

import edu.uta.futureye.core.DOF;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.TriAreaCoord;
import edu.uta.futureye.core.Vertex;
import edu.uta.futureye.core.intf.FiniteElement;
import edu.uta.futureye.function.basic.FX;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.util.container.VertexList;

public class FELinearTriangle implements FiniteElement {
	TriAreaCoord coord;
	
	//Construct a function with the coordinate of points in an element as parameters
	String[] argsOrder;// = new String[]{"x1","x2","x3","y1","y2","y3","r","s","t"};
	
	public int nDOFs = 3;
	MathFunc[] shapeFuncs = new MathFunc[nDOFs];

	public FELinearTriangle() {
		FX x1 = new FX("x1");
		FX x2 = new FX("x2");
		FX x3 = new FX("x3");
		FX y1 = new FX("y1");
		FX y2 = new FX("y2");
		FX y3 = new FX("y3");
		
		coord = new TriAreaCoord(x1, x2, x3, y1, y2, y3);
		MathFunc r = coord.getCoordR();
		MathFunc s = coord.getCoordS();
		
		argsOrder = new String[]{x1, x2, x3, y1, y2, y3, r, s, "t"};
		
		//shape functions
		shapeFuncs[0] = r;
		shapeFuncs[1] = s;
		shapeFuncs[2] = 1 - r - s;
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
		return this.coord.Jacobian();
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
