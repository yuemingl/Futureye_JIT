package edu.uta.futureye.lib.element;

import edu.uta.futureye.core.DOF;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Vertex;
import edu.uta.futureye.lib.shapefun.SFBilinearLocal2D;
import edu.uta.futureye.util.container.VertexList;

public class FEBilinearRectangle implements FiniteElementType {
	SFBilinearLocal2D[] shapeFun = new SFBilinearLocal2D[4];

	public FEBilinearRectangle() {
		shapeFun[0] = new SFBilinearLocal2D(1);
		shapeFun[1] = new SFBilinearLocal2D(2);
		shapeFun[2] = new SFBilinearLocal2D(3);
		shapeFun[3] = new SFBilinearLocal2D(4);
	}
	
	/**
	 * Assign degree of freedom to element
	 * @param e
	 */
	@Override
	public void assignTo(Element e) {
		//Asign degree of freedom to element
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
		return 4;
	}

	@Override
	public int getVectorShapeFunctionDim() {
		throw new UnsupportedOperationException();
	}

}
