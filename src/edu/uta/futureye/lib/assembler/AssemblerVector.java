package edu.uta.futureye.lib.assembler;

import edu.uta.futureye.algebra.SparseBlockMatrix;
import edu.uta.futureye.algebra.SparseBlockVector;
import edu.uta.futureye.algebra.SparseMatrix;
import edu.uta.futureye.algebra.SparseVector;
import edu.uta.futureye.algebra.intf.BlockMatrix;
import edu.uta.futureye.algebra.intf.BlockVector;
import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.core.DOF;
import edu.uta.futureye.core.DOFOrder;
import edu.uta.futureye.core.EdgeLocal;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.FaceLocal;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.core.NodeType;
import edu.uta.futureye.core.geometry.GeoEntity;
import edu.uta.futureye.core.intf.Assembler;
import edu.uta.futureye.core.intf.WeakForm;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.function.intf.VectorFunction;
import edu.uta.futureye.function.intf.VectorShapeFunction;
import edu.uta.futureye.lib.element.FiniteElementType;
import edu.uta.futureye.util.container.DOFList;
import edu.uta.futureye.util.container.ElementList;
import edu.uta.futureye.util.container.NodeList;

public class AssemblerVector implements Assembler {
	protected Mesh mesh;
	protected WeakForm weakForm;
	protected BlockMatrix globalStiff;
	protected BlockVector globalLoad;
	
	public AssemblerVector(Mesh mesh, WeakForm weakForm, FiniteElementType feType) {
		this.mesh = mesh;
		this.weakForm = weakForm;
		
		ElementList eList = mesh.getElementList();
		int nEle = eList.size();
		NodeList nodes = mesh.getNodeList();
		int nNode = nodes.size();
		
//		Element e = eList.at(1);
//		DOFList DOFs = e.getAllDOFList(DOFOrder.NEFV);
//		VectorFunction vsf = DOFs.at(1).getVSF();
//		int blockDim = vsf.getDim();
		int blockDim = feType.getVectorShapeFunctionDim();
		
		//TODO
		int[] dims = {0,nNode,nNode,mesh.nVertex};
		globalStiff = new SparseBlockMatrix(blockDim,blockDim);
		globalLoad = new SparseBlockVector(blockDim);
		for(int i=1;i<=blockDim;i++) {
			for(int j=1;j<=blockDim;j++) {
				globalStiff.setBlock(i, j, 
						new SparseMatrix(dims[i],dims[j]));
			}
			globalLoad.setBlock(i, 
					new SparseVector(dims[i]));
		}
		
	}
	
	@Override
	public void assemble() {
		ElementList eList = mesh.getElementList();
		int nEle = eList.size();
		for(int i=1; i<=nEle; i++) {
			eList.at(i).adjustVerticeToCounterClockwise();
			assembleGlobal(eList.at(i),	globalStiff,globalLoad);
			if(i%3000==0)
				System.out.println("Assemble..."+
						String.format("%.0f%%", 100.0*i/nEle));
		}
		//procHangingNode(mesh);
	}
	

	/**
	 * 从单元e合成全局矩阵和向量
	 * @param e
	 * @param stiff
	 * @param load
	 */
	public void assembleGlobal(Element e, Matrix stiff, Vector load) {
		DOFList DOFs = e.getAllDOFList(DOFOrder.NEFV);
		int nDOFs = DOFs.size();
		
		//Update Jacobin on e
		e.updateJacobin();
		
		//形函数计算需要和单元关联
		for(int i=1;i<=nDOFs;i++) {
			DOFs.at(i).getVSF().asignElement(e);
		}
		
		//所有自由度双循环
		for(int i=1;i<=nDOFs;i++) {
			DOF dofI = DOFs.at(i);
			VectorShapeFunction sfI = dofI.getVSF();
			int nLocalRow = dofI.getLocalIndex();
			int nGlobalRow = dofI.getGlobalIndex();
			int vvfIndexI = dofI.getVvfIndex();
			for(int j=1;j<=nDOFs;j++) {
				DOF dofJ = DOFs.at(j);
				VectorShapeFunction sfJ = dofJ.getVSF();
				int nLocalCol = dofJ.getLocalIndex();
				int nGlobalCol = dofJ.getGlobalIndex();
				int vvfIndexJ = dofJ.getVvfIndex();
				//???加上有错误，为什么？
				//if(vvfIndexI == vvfIndexJ) { //不相等合成结果是0，不用计算，e.g. Stokes:(u v p)
					//Local stiff matrix
					//注意顺序，内循环test函数不变，trial函数循环
					weakForm.setShapeFunction(sfJ,nLocalCol, sfI,nLocalRow); 
					Function lhs = weakForm.leftHandSide(e, WeakForm.ItemType.Domain);
					double lhsVal = weakForm.integrate(e, lhs);
					stiff.add(nGlobalRow, nGlobalCol, lhsVal);
				//}
			}
			//Local load vector
			weakForm.setShapeFunction(null,0,sfI,nLocalRow);
			Function rhs = weakForm.rightHandSide(e, WeakForm.ItemType.Domain);
			double rhsVal = weakForm.integrate(e, rhs);
			load.add(nGlobalRow, rhsVal);
		}
		
		if(e.isBorderElement()) {
			ElementList beList = e.getBorderElements();
			for(int n=1;n<=beList.size();n++) {
				Element be = beList.at(n);
				be.updateJacobin();
				DOFList beDOFs = be.getAllDOFList(DOFOrder.NEFV);
				int nBeDOF = beDOFs.size();
				
				//形函数计算需要和单元关联
				for(int i=1;i<=nBeDOF;i++) {
					beDOFs.at(i).getVSF().asignElement(e);
				}
				
				//所有自由度双循环
				for(int i=1;i<=nBeDOF;i++) {
					DOF dofI = beDOFs.at(i);
					VectorShapeFunction sfI = dofI.getVSF();
					int nLocalRow = dofI.getLocalIndex();
					int nGlobalRow = dofI.getGlobalIndex();
					int vvfIndexI = dofI.getVvfIndex();
					for(int j=1;j<=nBeDOF;j++) {
						DOF dofJ = beDOFs.at(j);
						VectorShapeFunction sfJ = dofJ.getVSF();
						int nLocalCol = dofJ.getLocalIndex();
						int nGlobalCol = dofJ.getGlobalIndex();
						int vvfIndexJ = dofJ.getVvfIndex();
						
						if(vvfIndexI == vvfIndexJ) {
							//Check node type
							NodeType nodeType = be.getBorderNodeType(vvfIndexI);
							if(nodeType == NodeType.Neumann || nodeType == NodeType.Robin) {
								//Local stiff matrix for border
								//注意顺序，内循环test函数不变，trial函数循环
								weakForm.setShapeFunction(sfJ,nLocalCol, sfI,nLocalRow);
								Function lhsBr = weakForm.leftHandSide(be, WeakForm.ItemType.Border);
								double lhsBrVal = weakForm.integrate(be, lhsBr);
								stiff.add(nGlobalRow, nGlobalCol, lhsBrVal);
							}
						}
					}
					//Local load vector for border
					weakForm.setShapeFunction(null,0,sfI,nLocalRow);
					Function rhsBr = weakForm.rightHandSide(be, WeakForm.ItemType.Border);
					double rhsBrVal = weakForm.integrate(be, rhsBr);
					load.add(nGlobalRow, rhsBrVal);
				}
			}
		}
	}
	
	@Override
	public Vector getLoadVector() {
		return this.globalLoad;
	}

	@Override
	public Matrix getStiffnessMatrix() {
		return this.globalStiff;
	}

	@Override
	public void imposeDirichletCondition(Function diri) {
		throw new UnsupportedOperationException();

	}
	
	
	protected void setDirichlet(int matIndex, double value) {
		int row = matIndex;
		int col = matIndex;
		this.globalStiff.set(row, col, 1.0);
		this.globalLoad.set(row,value);
		for(int r=1;r<=this.globalStiff.getRowDim();r++) {
			if(r != row) {
				this.globalLoad.add(r,-this.globalStiff.get(r, col)*value);
				this.globalStiff.set(r, col, 0.0);
			}
		}
		for(int c=1;c<=this.globalStiff.getColDim();c++) {
			if(c != col) {
				this.globalStiff.set(row, c, 0.0);
			}
		}
	}
	
	/**
	 * 向量值问题的Dirichlet条件在整个矩阵上施加，而不是分块矩阵上
	 * 
	 */
	@Override
	public void imposeDirichletCondition(VectorFunction diri) {
		ElementList eList = mesh.getElementList();
		for(int i=1;i<=eList.size();i++) {
			Element e = eList.at(i);
			DOFList DOFs = e.getAllDOFList(DOFOrder.NEFV);
			for(int j=1;j<=DOFs.size();j++) {
				DOF dof = DOFs.at(j);
				GeoEntity ge = dof.getOwner();
				int vvfIndex = dof.getVvfIndex();
				Function fdiri = diri.get(vvfIndex);
				if(ge instanceof Node) {
					Node n = (Node)ge;
					if(n.getNodeType(vvfIndex) == NodeType.Dirichlet) {
						Variable v = Variable.createFrom(fdiri, n, 0);
						setDirichlet(dof.getGlobalIndex(),fdiri.value(v));
					}
				} else if(ge instanceof EdgeLocal) {
					EdgeLocal edge = (EdgeLocal)ge;
					if(edge.getBorderType(vvfIndex) == NodeType.Dirichlet) {
						//TODO 以边的那个顶点取值？中点？
						//Variable v = Variable.createFrom(fdiri, ?, 0);
					}
					
				} else if(ge instanceof FaceLocal) {
					FaceLocal face = (FaceLocal)ge;
					if(face.getBorderType(vvfIndex) == NodeType.Dirichlet) {
						//TODO
					}
				}
			}
			
		}
	}		
}
