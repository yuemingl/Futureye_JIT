package edu.uta.futureye.lib.element;

import edu.uta.futureye.core.DOF;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.lib.shapefun.SFQuadraticLocal2DFast;

public class FEQuadraticTriangle implements FiniteElementType {
	protected SFQuadraticLocal2DFast[] shapeFun = new SFQuadraticLocal2DFast[6];
//	protected SFQuadraticLocal2D[] shapeFun = new SFQuadraticLocal2D[6];

	public FEQuadraticTriangle() {
		for(int i=1;i<=6;i++)
			shapeFun[i-1] = new SFQuadraticLocal2DFast(i);
//		for(int i=1;i<=6;i++)
//			shapeFun[i-1] = new SFQuadraticLocal2D(i);
	}
	
	@Override
	public void assignTo(Element e) {
		//Assign shape function to DOF
		for(int j=1;j<=e.nodes.size();j++) {
			//Asign shape function to DOF
			DOF dof = new DOF(
					j,//Local DOF index
					e.nodes.at(j).globalIndex,//Global DOF index, take global node index
					shapeFun[j-1]//Shape function 
					         );
			e.addNodeDOF(j, dof);
		}
	}

	@Override
	public int getDOFNumOnElement(int vsfDim) {
		return 6;
	}

	@Override
	public int getVectorShapeFunctionDim() {
		throw new UnsupportedOperationException();
	}
}
