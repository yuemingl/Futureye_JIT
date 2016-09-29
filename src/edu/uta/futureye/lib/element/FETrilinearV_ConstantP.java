package edu.uta.futureye.lib.element;

import edu.uta.futureye.core.DOF;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.lib.shapefun.TrilinearV_ConstantP;
import edu.uta.futureye.util.FutureyeException;

public class FETrilinearV_ConstantP implements FiniteElementType {
	protected static TrilinearV_ConstantP[] shapeFun = new TrilinearV_ConstantP[25];
	protected int nTotalNodes = -1;
	//p自由度计数器
	protected int nDOF_p = -1;
	
	public int getVectorShapeFunctionDim() {
		return 4;
	}
	
	public int getDOFNumOnElement(int vsfDim) {
		if(vsfDim <= 3)
			return 8;
		else
			return 1;
	}
	
	public FETrilinearV_ConstantP() {
		for(int i=0;i<25;i++)
			shapeFun[i] = new TrilinearV_ConstantP(i+1);
	}
	
	/**
	 * Assign degree of freedom to element
	 * @param e
	 */
	public void assignTo(Element e) {
		if(nTotalNodes == -1 || nDOF_p == -1) {
			FutureyeException ex = new FutureyeException("Call initDOFIndex() first!");
			ex.printStackTrace();
			System.exit(-1);
		}
		//单元结点数
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
			dof_u1.setVVFComponent(1);
			DOF dof_u2 = new DOF(
					nNode+j,//Local DOF index
					//Global DOF index, take this.nTotalNodes + global node index
					this.nTotalNodes+e.nodes.at(j).globalIndex,
					shapeFun[nNode+j-1]//Shape function 
					         );
			dof_u2.setVVFComponent(2);
			DOF dof_u3 = new DOF(
					nNode*2+j,//Local DOF index
					//Global DOF index, take this.nTotalNodes*2 + global node index
					this.nTotalNodes*2+e.nodes.at(j).globalIndex,
					shapeFun[nNode*2+j-1]//Shape function 
					         );
			dof_u3.setVVFComponent(3);
			e.addNodeDOF(j, dof_u1);
			e.addNodeDOF(j, dof_u2);
			e.addNodeDOF(j, dof_u3);
		}
		
		//Assign shape function to DOF
		DOF dof = new DOF(
					nNode*3+1, //Local DOF index
					this.nTotalNodes*3+this.nDOF_p, //Global DOF index for Pressure
					shapeFun[nNode*3] //Shape function 
					);
		this.nDOF_p++;
		dof.setVVFComponent(4);	
		e.addVolumeDOF(dof);
	}
	
	@Override
	public int getDOFNumOnMesh(Mesh mesh, int vsfDim) {
		if(vsfDim<=3)
			return mesh.getNodeList().size();
		else
			return mesh.getElementList().size();
	}

	@Override
	public void initDOFIndexGenerator(Mesh mesh) {
		this.nTotalNodes = mesh.getNodeList().size();
		nDOF_p = 1;
	}	
}
