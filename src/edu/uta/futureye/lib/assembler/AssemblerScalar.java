package edu.uta.futureye.lib.assembler;

import java.util.Map;
import java.util.Map.Entry;

import symjava.bytecode.BytecodeFunc;
import symjava.symbolic.Expr;
import symjava.symbolic.Func;
import edu.uta.futureye.algebra.SparseMatrixRowMajor;
import edu.uta.futureye.algebra.SparseVectorHashMap;
import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.SparseMatrix;
import edu.uta.futureye.algebra.intf.SparseVector;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.core.DOF;
import edu.uta.futureye.core.DOFOrder;
import edu.uta.futureye.core.Edge;
import edu.uta.futureye.core.EdgeLocal;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Face;
import edu.uta.futureye.core.FaceLocal;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.core.NodeRefined;
import edu.uta.futureye.core.NodeType;
import edu.uta.futureye.core.Volume;
import edu.uta.futureye.core.geometry.GeoEntity;
import edu.uta.futureye.core.geometry.GeoEntity2D;
import edu.uta.futureye.core.geometry.GeoEntity3D;
import edu.uta.futureye.core.intf.Assembler;
import edu.uta.futureye.core.intf.WeakForm;
import edu.uta.futureye.function.FEMFunc;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.intf.MathFun;
import edu.uta.futureye.function.intf.VectorFunction;
import edu.uta.futureye.util.FutureyeException;
import edu.uta.futureye.util.container.DOFList;
import edu.uta.futureye.util.container.ElementList;
import edu.uta.futureye.util.container.VertexList;
import edu.uta.futureye.core.intf.WeakForm.ItemType;

public class AssemblerScalar implements Assembler {
	private int status = 0;
	protected Mesh mesh;
	protected WeakForm weakForm;
	protected SparseMatrix globalStiff;
	protected SparseVector globalLoad;
	private boolean printInfo = true;

	public AssemblerScalar(Mesh mesh, WeakForm weakForm) {
		this.mesh = mesh;
		this.weakForm = weakForm;
		
		int dim = mesh.getNodeList().size();
		globalStiff = new SparseMatrixRowMajor(dim,dim);
		globalLoad = new SparseVectorHashMap(dim);

	}
	
	@Override
	public SparseMatrix getStiffnessMatrix() {
		if(status == 0)
			throw new FutureyeException("Call assemble() function first!");
		return this.globalStiff;
	}
	
	@Override
	public SparseVector getLoadVector() {
		if(status == 0)
			throw new FutureyeException("Call assemble() function first!");
		return this.globalLoad;
	}
	
	@Override
	public void assemble() {
		assemble(true);
	}
	
	public void assemble(boolean procHangingNode) {
		status = 1;
		ElementList eList = mesh.getElementList();
		int nEle = eList.size();
		int nProgress = 20;
		if(printInfo) {
			System.out.print("Assemble[0%");
			for(int i=0;i<nProgress-6;i++)
				System.out.print("-");
			System.out.println("100%]");
			System.out.print("Progress[");
		}
		
		int nPS = nEle/nProgress;
		int nProgressPrint = 0;
		for(int i=1; i<=nEle; i++) {
			eList.at(i).adjustVerticeToCounterClockwise();
			//TODO
//			if(eList.at(i).adjustVerticeToCounterClockwise()) {
//				throw new FutureyeException("adjustVerticeToCounterClockwise");
//			}
			assembleGlobal(eList.at(i),	globalStiff,globalLoad);
			if(printInfo) {
				if(i%nPS==0) {
					nProgressPrint++;
					System.out.print("*");
				}
			}
		}
		if(printInfo) {
			if(nProgressPrint<nProgress)
				System.out.print("*");
			System.out.println("]Done!");
		}
		
		if(procHangingNode)
			procHangingNode(mesh);
	}
	
	public void printInfo(boolean flag) {
		this.printInfo = flag;
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
	
//	@Override
//	public void imposeDirichletCondition(Function diri) {
//		status = 2;
//		NodeList nList = mesh.getNodeList();
//		for(int i=1;i<=nList.size();i++) {
//			Node n = nList.at(i);
//			if(n.getNodeType() == NodeType.Dirichlet) {
//				Variable v = Variable.createFrom(diri, n, n.globalIndex);
//				this.globalStiff.set(n.globalIndex, n.globalIndex, 1.0);
//				double val = diri.value(v);
//				this.globalLoad.set(n.globalIndex, val);
//				for(int j=1;j<=this.globalStiff.getRowDim();j++) {
//					if(j!=n.globalIndex) {
//						//TODO 行列都需要置零
//						this.globalLoad.add(j, -this.globalStiff.get(j, n.globalIndex)*val);
//						this.globalStiff.set(j, n.globalIndex, 0.0);
//						this.globalStiff.set(n.globalIndex, j, 0.0);
//					}
//				}
//			}
//		}
//	}
	
	//2011/11/28 modified according to vector valued case
	@Override
	public void imposeDirichletCondition(Expr diri) {
		Func fdiri = new Func("diri", diri);
		BytecodeFunc bDiri = fdiri.toBytecodeFunc();
		ElementList eList = mesh.getElementList();
		for(int i=1;i<=eList.size();i++) {
			Element e = eList.at(i);
			DOFList DOFs = e.getAllDOFList(DOFOrder.NEFV);
			for(int j=1;j<=DOFs.size();j++) {
				DOF dof = DOFs.at(j);
				GeoEntity ge = dof.getOwner();
				if(ge instanceof Node) {
					Node n = (Node)ge;
					if(n.getNodeType() == NodeType.Dirichlet) {
						//Variable v = Variable.createFrom(diri, n, n.globalIndex); //bugfix 11/27/2013 Variable.createFrom(diri, n, 0);
						setDirichlet(dof.getGlobalIndex(), bDiri.apply(n.coords()));
					}
				} else if(ge instanceof EdgeLocal) {
					//2D单元（面）其中的局部边上的自由度
					EdgeLocal edge = (EdgeLocal)ge;
					if(edge.getBorderType() == NodeType.Dirichlet) {
						//TODO 以边的那个顶点取值？中点？
						//Variable v = Variable.createFrom(fdiri, ?, 0);
					}
					
				} else if(ge instanceof FaceLocal) {
					//3D单元（体）其中的局部面上的自由度
					FaceLocal face = (FaceLocal)ge;
					if(face.getBorderType() == NodeType.Dirichlet) {
						//TODO
					}
				} else if(ge instanceof Edge) {
					//1D单元（线段）上的自由度，其Dirichlet边界用结点来计算推出，而不需要专门标记单元
					VertexList vs = ((GeoEntity2D) ge).getVertices();
					for(int k=1;k<=vs.size();k++) {
						Node n = vs.at(k).globalNode();
						if(NodeType.Dirichlet == n.getNodeType()) {
							//Variable v = Variable.createFrom(diri, n, 0);
							//setDirichlet(dof.getGlobalIndex(),diri.apply(v));
							setDirichlet(dof.getGlobalIndex(), bDiri.apply(n.coords()));
						}
					}
				} else if(ge instanceof Face) {
					//2D单元（面）上的自由度，其Dirichlet边界用结点来计算推出，而不需要专门标记单元
					
					VertexList vs = ((GeoEntity2D) ge).getVertices();
					for(int k=1;k<=vs.size();k++) {
						Node n = vs.at(k).globalNode();
						if(NodeType.Dirichlet == n.getNodeType()) {
							//Variable v = Variable.createFrom(diri, n, 0);
							//setDirichlet(dof.getGlobalIndex(),diri.apply(v));
							setDirichlet(dof.getGlobalIndex(), bDiri.apply(n.coords()));
						}
					}
				} else if(ge instanceof Volume) {
					//3D单元（体）上的自由度，其Dirichlet边界用结点来计算推出，而不需要专门标记单元
					VertexList vs = ((GeoEntity3D) ge).getVertices();
					for(int k=1;k<=vs.size();k++) {
						Node n = vs.at(k).globalNode();
						if(NodeType.Dirichlet == n.getNodeType()) {
							//Variable v = Variable.createFrom(diri, n, 0);
							//setDirichlet(dof.getGlobalIndex(),diri.apply(v));
							setDirichlet(dof.getGlobalIndex(), bDiri.apply(n.coords()));
						}
					}
				}
			}
		}
	}
	
	/**
	 * 从单元e合成全局矩阵和向量
	 * @param e
	 * @param stiff
	 * @param load
	 */
	public void assembleGlobal(Element e, Matrix stiff, Vector load) {
		status = 3;
		DOFList DOFs = e.getAllDOFList(DOFOrder.NEFV);
		int nDOFs = DOFs.size();
		
		//Update Jacobin on e
		e.updateJacobin();
		//System.out.println(e.getJacobin().value(new Variable("r",0).set("s", 0)));
		
		//形函数计算需要和单元关联
		for(int i=1;i<=nDOFs;i++) {
			DOFs.at(i).getSSF().assignElement(e);
		}
		
		weakForm.preProcess(e);
		
		//所有自由度双循环
		for(int i=1;i<=nDOFs;i++) {
			DOF dofI = DOFs.at(i);
			int nGlobalRow = dofI.getGlobalIndex();
			for(int j=1;j<=nDOFs;j++) {
				DOF dofJ = DOFs.at(j);
				int nGlobalCol = dofJ.getGlobalIndex();
				//Local stiff matrix
				//注意顺序，内循环test基函数不变，trial基函数循环
				weakForm.setDOF(dofJ, dofI); 
				
//				MathFun lhs = weakForm.leftHandSide(e, ItemType.Domain);
//				double lhsVal = weakForm.integrate(e, lhs);
//				stiff.add(nGlobalRow, nGlobalCol, lhsVal);
				
				FEMFunc lhs = weakForm.leftHandSide(e, ItemType.Domain);
				stiff.add(nGlobalRow, nGlobalCol, lhs.apply(e.getNodeCoords()));
				
				
			}
			//Local load vector
			weakForm.setDOF(null,dofI);
//			MathFun rhs = weakForm.rightHandSide(e, ItemType.Domain);
//			double rhsVal = weakForm.integrate(e, rhs);
//			load.add(nGlobalRow, rhsVal);
			FEMFunc rhs = weakForm.rightHandSide(e, ItemType.Domain);
			load.add(nGlobalRow, rhs.apply(e.getNodeCoords()));

		}
		
		if(e.isBorderElement()) {
			ElementList beList = e.getBorderElements();
			for(int n=1;n<=beList.size();n++) {
				Element be = beList.at(n);
				
				//Check node type
				NodeType nodeType = be.getBorderNodeType();
				if(nodeType == NodeType.Neumann || nodeType == NodeType.Robin) {
					be.updateJacobin();
					weakForm.preProcess(be);
					DOFList beDOFs = be.getAllDOFList(DOFOrder.NEFV);
					int nBeDOF = beDOFs.size();
					
					//形函数计算需要和单元关联
					for(int i=1;i<=nBeDOF;i++) {
						beDOFs.at(i).getSSF().assignElement(be);
					}
					
					//所有自由度双循环
					for(int i=1;i<=nBeDOF;i++) {
						DOF dofI = beDOFs.at(i);
						int nGlobalRow = dofI.getGlobalIndex();
						for(int j=1;j<=nBeDOF;j++) {
							DOF dofJ = beDOFs.at(j);
							int nGlobalCol = dofJ.getGlobalIndex();
							//Local stiff matrix for border
							//注意顺序，内循环test函数不变，trial函数循环
							weakForm.setDOF(dofJ, dofI);
							FEMFunc lhsBr = weakForm.leftHandSide(be, ItemType.Border);
							stiff.add(nGlobalRow, nGlobalCol, lhsBr.apply(be.getNodeCoords()));
						}
						//Local load vector for border
						weakForm.setDOF(null, dofI);
						FEMFunc rhsBr = weakForm.rightHandSide(be, ItemType.Border);
						load.add(nGlobalRow, rhsBr.apply(be.getNodeCoords()));
					}
				}
			}
		}
	}
	
	/**
  	 * 从单元e合成局部矩阵和向量
  	 * @param e
	 * @param localStiff (I/O): 
	 *   local stiff matrix corresponds to the integration part
	 *   in element e
	 * @param localStiffBorder (I/O): 
	 *   local stiff matrix corresponds to the integration part 
	 *   on the border of element e
	 * @param localLoad (I/O): 
	 *   local load vector
	 * @param localLoadBorder (I/O): 
	 *   local load vector for border
	 */
	public void assembleLocal(Element e, 
			SparseMatrix localStiff, SparseMatrix localStiffBorder, 
			SparseVector localLoad, SparseVector localLoadBorder) {
		status = 4;
		
		localStiff.clearAll();
		localLoad.clearAll();
		
		DOFList DOFs = e.getAllDOFList(DOFOrder.NEFV);
		int nDOFs = DOFs.size();
		
		localStiff.setRowDim(nDOFs);
		localStiff.setColDim(nDOFs);
		localLoad.setDim(nDOFs);
		
		//Update Jacobin on e
		e.updateJacobin();
		
		//形函数计算需要和单元关联
		for(int i=1;i<=nDOFs;i++) {
			DOFs.at(i).getSSF().assignElement(e);
		}
		
		//所有自由度双循环
		for(int i=1;i<=nDOFs;i++) {
			DOF dofI = DOFs.at(i);
			int nLocalRow = dofI.getLocalIndex();
			for(int j=1;j<=nDOFs;j++) {
				DOF dofJ = DOFs.at(j);
				int nLocalCol = dofJ.getLocalIndex();
				//Local stiff matrix
				//注意顺序，内循环test函数不变，trial函数循环
				weakForm.setDOF(dofJ, dofI);
				FEMFunc lhs = weakForm.leftHandSide(e, WeakForm.ItemType.Domain);
				localStiff.add(nLocalRow, nLocalCol, lhs.apply(e.getNodeCoords()));
			}
			//Local load vector
			weakForm.setDOF(null, dofI);
			FEMFunc rhs = weakForm.rightHandSide(e, WeakForm.ItemType.Domain);
			localLoad.add(nLocalRow, rhs.apply(e.getNodeCoords()));
		}
		
		if(e.isBorderElement()) {
			ElementList beList = e.getBorderElements();
			for(int n=1;n<=beList.size();n++) {
				Element be = beList.at(n);
				//Check node type
				NodeType nodeType = be.getBorderNodeType();
				if(nodeType == NodeType.Neumann || nodeType == NodeType.Robin) {
					be.updateJacobin();
					DOFList beDOFs = be.getAllDOFList(DOFOrder.NEFV);
					int nBeDOF = beDOFs.size();
					
					localStiffBorder.setRowDim(nDOFs);
					localStiffBorder.setColDim(nDOFs);
					//形函数计算需要和单元关联
					for(int i=1;i<=nBeDOF;i++) {
						beDOFs.at(i).getSSF().assignElement(e);
					}
					
					//所有自由度双循环
					for(int i=1;i<=nBeDOF;i++) {
						DOF dofI = beDOFs.at(i);
						int nLocalRow = dofI.getLocalIndex();
						for(int j=1;j<=nBeDOF;j++) {
							DOF dofJ = beDOFs.at(j);
							int nLocalCol = dofJ.getLocalIndex();
							//Local stiff matrix for border
							//注意顺序，内循环test函数不变，trial函数循环
							weakForm.setDOF(dofJ, dofI);
							FEMFunc lhsBr = weakForm.leftHandSide(be, WeakForm.ItemType.Border);
							localStiffBorder.add(nLocalRow, nLocalCol, lhsBr.apply(be.getNodeCoords()));
						}
						//Local load vector for border
						weakForm.setDOF(null, dofI);
						FEMFunc rhsBr = weakForm.rightHandSide(be, WeakForm.ItemType.Border);
						localLoadBorder.add(nLocalRow, rhsBr.apply(be.getNodeCoords()));
					}
				}
			}
		}
	}
	
	//二维：刚度矩阵增加hanging node约束系数
	// nh - 0.5*n1 - 0.5*n2 = 0
	//不能去掉， hanging node 自己在刚度矩阵合成的时候没有系数，因此需要为自己设置系数
	public void procHangingNode(Mesh mesh) {
		status = 5;
		for(int i=1;i<=mesh.getNodeList().size();i++) {
			Node node = mesh.getNodeList().at(i);
			if(node instanceof NodeRefined) {
				NodeRefined nRefined = (NodeRefined)node;
				if(nRefined.isHangingNode()) {
					globalStiff.set(nRefined.globalIndex, nRefined.globalIndex, 1.0);
					globalStiff.set(nRefined.globalIndex,
							nRefined.constrainNodes.at(1).globalIndex,-0.5);
					globalStiff.set(nRefined.globalIndex,
							nRefined.constrainNodes.at(2).globalIndex,-0.5);
				}
			}
		}
	}

	public void plusToGlobalLoad(Element e, Vector local) {
		for(int i=1; i<=local.getDim(); i++) {
			int ngRow = e.local2GlobalDOFIndex(i);
			this.globalLoad.add(ngRow, local.get(i));
		}
	}
	
	public void plusToGlobalStriff(Element e, SparseMatrix local) {
		Map<Integer,Map<Integer,Double>> m = local.getAll();
		for(Entry<Integer,Map<Integer,Double>> er : m.entrySet()) {
			int nRow = er.getKey();
			int ngRow = e.local2GlobalDOFIndex(nRow);
			Map<Integer,Double> row = er.getValue();
			for(Entry<Integer,Double> ec : row.entrySet()) {
				int nCol = ec.getKey();
				double val = ec.getValue();
				int ngCol = e.local2GlobalDOFIndex(nCol);
//				System.out.println(
//						"Local("+nRow+","+nCol+") -> Global("+ngRow+","+ngCol+")"
//						);
				this.globalStiff.add(ngRow, ngCol, val);
			}
		}
	}
	
	@Override
	public void imposeDirichletCondition(VectorFunction diri) {
		throw new UnsupportedOperationException();
	}	
}
