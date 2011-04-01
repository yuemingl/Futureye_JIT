package edu.uta.futureye.lib.element;

import edu.uta.futureye.core.DOF;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.lib.shapefun.QuadraticV_LinearP;
import edu.uta.futureye.util.FutureyeException;
import edu.uta.futureye.util.container.VertexList;

public class FEQuadraticV_LinearP implements FiniteElementType {
	protected static QuadraticV_LinearP[] shapeFun = new QuadraticV_LinearP[15];
	protected int nTotalNodes = -1;
	protected int nDOF_p = -1;
	
	public int getVectorShapeFunctionDim() {
		return 3;
	}
	
	public int getDOFNumOnElement(int vsfDim) {
		if(vsfDim == 1 || vsfDim == 2)
			return 6;
		else
			return 3;
	}
	
	
	public FEQuadraticV_LinearP() {
		for(int i=0;i<15;i++)
			shapeFun[i] = new QuadraticV_LinearP(i+1);
	}
	
	public void initDOFIndexGenerator(int nTotalNodes) {
		this.nTotalNodes = nTotalNodes;
		nDOF_p = 1;
	}
	
	/**
	 * Assign degree of freedom to element
	 * @param e
	 */
	public void assign(Element e) {
		if(nTotalNodes == -1 || nDOF_p == -1) {
			FutureyeException ex = new FutureyeException("Call initDOFIndex() first!");
			ex.printStackTrace();
			System.exit(-1);
		}
		int nNode = e.nodes.size();
		//Assign shape function to DOF
		for(int j=1;j<=nNode;j++) {
			//Asign shape function to DOF
			DOF dof_u1 = new DOF(
					j,//Local DOF index
					//Global DOF index, take global node index
					e.nodes.at(j).globalIndex,
					shapeFun[j-1]//Shape function 
					         );
			dof_u1.setVvfIndex(1);
			DOF dof_u2 = new DOF(
					nNode+j,//Local DOF index
					//Global DOF index, take this.nTotalNodes + global node index
					this.nTotalNodes+e.nodes.at(j).globalIndex,
					shapeFun[nNode+j-1]//Shape function 
					         );
			dof_u2.setVvfIndex(2);
			e.addNodeDOF(j, dof_u1);
			e.addNodeDOF(j, dof_u2);
		}
		VertexList vertices = e.vertices();
		for(int j=1;j<=vertices.size();j++) {
			//Assign shape function to DOF
			DOF dof = new DOF(
						2*nNode+j, //Local DOF index
						//this.nTotalNodes*2+nDOF_p, //Global DOF index for Pressure
						this.nTotalNodes*2+vertices.at(j).globalNode().globalIndex, //Global DOF index for Pressure
						shapeFun[2*nNode+j-1] //Shape function 
						);
			dof.setVvfIndex(3);
			//System.out.println(this.nTotalNodes*2+nDOF_p);
			nDOF_p++;
			e.addNodeDOF(j, dof);
		}
	}
}
