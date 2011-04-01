package edu.uta.futureye.lib.element;

import edu.uta.futureye.core.DOF;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Vertex;
import edu.uta.futureye.lib.shapefun.SFLinearLocal2D;
import edu.uta.futureye.util.container.VertexList;

/**
 * Linear Triangle element for 2D
 * 2维三角形线性单元
 * 
 * @author liuyueming
 *
 */
public class FELinearTriangle implements FiniteElementType {
	protected static SFLinearLocal2D[] shapeFun = new SFLinearLocal2D[3];
	
	public FELinearTriangle() {
		shapeFun[0] = new SFLinearLocal2D(1);
		shapeFun[1] = new SFLinearLocal2D(2);
		shapeFun[2] = new SFLinearLocal2D(3);
	}
	
	/**
	 * Assign degree of freedom to element
	 * @param e
	 */
	public void assign(Element e) {
		VertexList vertices = e.vertices();
		for(int j=1;j<=vertices.size();j++) {
			Vertex v = vertices.at(j);
			//Assign shape function to DOF
			DOF dof = new DOF(
						j, //Local DOF index
						v.globalNode().getIndex(), //Global DOF index, take global node index
						shapeFun[j-1] //Shape function 
						);
			e.addNodeDOF(j, dof);
		}
	}

	@Override
	public int getDOFNumOnElement(int vsfDim) {
		return 3;
	}

	@Override
	public int getVectorShapeFunctionDim() {
		throw new UnsupportedOperationException();
	}
	
}
