package edu.uta.futureye.lib.assembler;

import java.util.HashSet;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import edu.uta.futureye.algebra.FullMatrix;
import edu.uta.futureye.algebra.SparseBlockMatrix;
import edu.uta.futureye.algebra.SparseBlockVector;
import edu.uta.futureye.algebra.SparseMatrixRowMajor;
import edu.uta.futureye.algebra.SparseVectorHashMap;
import edu.uta.futureye.algebra.intf.Matrix;
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
import edu.uta.futureye.core.NodeType;
import edu.uta.futureye.core.Volume;
import edu.uta.futureye.core.geometry.GeoEntity;
import edu.uta.futureye.core.geometry.GeoEntity2D;
import edu.uta.futureye.core.geometry.GeoEntity3D;
import edu.uta.futureye.core.intf.AssemblerOld;
import edu.uta.futureye.core.intf.WeakFormOld;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.function.intf.VecMathFunc;
import edu.uta.futureye.lib.element.FiniteElementType;
import edu.uta.futureye.util.Constant;
import edu.uta.futureye.util.container.DOFList;
import edu.uta.futureye.util.container.ElementList;
import edu.uta.futureye.util.container.VertexList;

public class AssemblerVector implements AssemblerOld {
	protected Mesh mesh;
	protected WeakFormOld weakForm;
	protected SparseBlockMatrix globalStiff;
	protected SparseBlockVector globalLoad;
	
	/**
	 * 构造一个整体合成器
	 * 
	 * @param mesh 网格
	 * @param weakForm 弱形式
	 * @param feType 有限元类型
	 */
	public AssemblerVector(Mesh mesh, WeakFormOld weakForm, 
			FiniteElementType feType) {
		this.mesh = mesh;
		this.weakForm = weakForm;
		
		int blockDim = feType.getVectorShapeFunctionDim();
		//获取网格自由度总数
		int[] dims = new int[blockDim];
		for(int i=1;i<=blockDim;i++) {
			dims[i-1] = feType.getDOFNumOnMesh(mesh, i);
		}
		globalStiff = new SparseBlockMatrix(blockDim,blockDim);
		globalLoad = new SparseBlockVector(blockDim);
		for(int i=1;i<=blockDim;i++) {
			for(int j=1;j<=blockDim;j++) {
				globalStiff.setBlock(i, j, 
						new SparseMatrixRowMajor(dims[i-1],dims[j-1]));
			}
			globalLoad.setBlock(i, 
					new SparseVectorHashMap(dims[i-1]));
		}
		
	}
	
	@Override
	public void assemble() {
		ElementList eList = mesh.getElementList();
		int nEle = eList.size();
		int nProgress = 20;
		System.out.print("Assemble[0%");
		for(int i=0;i<nProgress-6;i++)
			System.out.print("-");
		System.out.println("100%]");
		
		System.out.print("Progress[");
		int nPS = nEle/nProgress;
		int nProgressPrint = 0;
		for(int i=1; i<=nEle; i++) {
			eList.at(i).adjustVerticeToCounterClockwise();
			assembleGlobal(eList.at(i),	globalStiff,globalLoad);
			if(i%nPS==0) {
				nProgressPrint++;
				System.out.print("*");
			}
		}
		if(nProgressPrint<nProgress)
			System.out.print("*");
		System.out.println("]Done!");
		
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
			DOFs.at(i).getVSF().assignElement(e);
		}
		
		weakForm.preProcess(e);

		//所有自由度双循环
		for(int i=1;i<=nDOFs;i++) {
			DOF dofI = DOFs.at(i);
			int nGlobalRow = dofI.getGlobalIndex();
			int nVVFCmptI = dofI.getVVFComponent();
			for(int j=1;j<=nDOFs;j++) {
				DOF dofJ = DOFs.at(j);
				int nGlobalCol = dofJ.getGlobalIndex();
				int nVVFCmptJ = dofJ.getVVFComponent();
				//不耦合的函数合成结果是0，不用计算，e.g. Stokes:(u v p)
				if(weakForm.isVVFComponentCoupled(nVVFCmptI,nVVFCmptJ)) { 
					//Local stiff matrix
					//注意顺序，内循环test函数不变，trial函数循环
					weakForm.setDOF(dofJ, dofI);
					MathFunc lhs = weakForm.leftHandSide(e, WeakFormOld.ItemType.Domain);
					double lhsVal = weakForm.integrate(e, lhs);
					stiff.add(nGlobalRow, nGlobalCol, lhsVal);
					//System.out.println(nVVFCmptI+"  "+nVVFCmptJ+"   "+lhsVal);
				}
			}
			//Local load vector
			weakForm.setDOF(null, dofI);
			MathFunc rhs = weakForm.rightHandSide(e, WeakFormOld.ItemType.Domain);
			double rhsVal = weakForm.integrate(e, rhs);
			load.add(nGlobalRow, rhsVal);
		}
		
		if(e.isBorderElement()) {
			ElementList beList = e.getBorderElements();
			for(int n=1;n<=beList.size();n++) {
				Element be = beList.at(n);
				
				be.updateJacobin();
				weakForm.preProcess(be);
				
				DOFList beDOFs = be.getAllDOFList(DOFOrder.NEFV);
				int nBeDOF = beDOFs.size();
				
				//形函数计算需要和单元关联
				for(int i=1;i<=nBeDOF;i++) {
					beDOFs.at(i).getVSF().assignElement(be);
				}
				
				//所有自由度双循环
				for(int i=1;i<=nBeDOF;i++) {
					DOF dofI = beDOFs.at(i);
					int nGlobalRow = dofI.getGlobalIndex();
					int nVVFCmptI = dofI.getVVFComponent();
					for(int j=1;j<=nBeDOF;j++) {
						DOF dofJ = beDOFs.at(j);
						int nGlobalCol = dofJ.getGlobalIndex();
						int nVVFCmptJ = dofJ.getVVFComponent();
						
						//不耦合的函数合成结果是0，不用计算，e.g. Stokes:(u v p)
						if(weakForm.isVVFComponentCoupled(nVVFCmptI,nVVFCmptJ)) { 
							//Check node type
							NodeType nodeType = be.getBorderNodeType(nVVFCmptI);
							if(nodeType == NodeType.Neumann || nodeType == NodeType.Robin) {
								//Local stiff matrix for border
								//注意顺序，内循环test函数不变，trial函数循环
								weakForm.setDOF(dofJ, dofI);
								MathFunc lhsBr = weakForm.leftHandSide(be, WeakFormOld.ItemType.Border);
								double lhsBrVal = weakForm.integrate(be, lhsBr);
								stiff.add(nGlobalRow, nGlobalCol, lhsBrVal);
							}
						}
					}
					//Check node type
					NodeType nodeType = be.getBorderNodeType(nVVFCmptI);
					if(nodeType == NodeType.Neumann || nodeType == NodeType.Robin) {
						//Local load vector for border
						weakForm.setDOF(null, dofI);
						MathFunc rhsBr = weakForm.rightHandSide(be, WeakFormOld.ItemType.Border);
						double rhsBrVal = weakForm.integrate(be, rhsBr);
						load.add(nGlobalRow, rhsBrVal);
					}
				}
			}
		}
	}
	
	@Override
	public SparseBlockVector getLoadVector() {
		return globalLoad;
	}

	@Override
	public SparseBlockMatrix getStiffnessMatrix() {
		return globalStiff;
	}

	@Override
	public void imposeDirichletCondition(MathFunc diri) {
		throw new UnsupportedOperationException();

	}
	
	
	protected void setDirichlet(FullMatrix stiff, int matIndex, double value) {
		int row = matIndex-1;
		int col = matIndex-1;
		double[][] data = stiff.getData();

		data[row][col] = 1.0;
		this.globalLoad.set(row+1,value);
		
		int nRow = stiff.getRowDim();
		int nCol = stiff.getColDim();
		if(Math.abs(value) > Matrix.zeroEps) {
			for(int r=0; r<nRow; r++) { 
				if(r != row) {
					//bugfix
					//this.globalLoad.add(r+1,-this.globalStiff.get(r+1, col+1)*value);
					this.globalLoad.add(r+1,-data[r][col]*value);
					data[r][col] = 0.0; //col列除对角元外，都置零
				}
			}
		} else {
			for(int r=0; r<nRow; r++) {
				if(r != row) {
					data[r][col] = 0.0; //col列除对角元外，都置零
				}
			}
		}
		for(int c=0; c<nCol; c++) {
			if(c != col) {
				data[row][c] = 0.0; //row行除对角元外，都置零
			}
		}
	}
	
	protected void setDirichlet(int matIndex, double value) {
		int row = matIndex;
		int col = matIndex;
		this.globalStiff.set(row, col, 1.0);
		this.globalLoad.set(row,value);
		int nRow = globalStiff.getRowDim();
		int nCol = globalStiff.getColDim();
		//value!=0
		if(Math.abs(value) > Matrix.zeroEps) {
			for(int r=1;r<=nRow;r++) {
				if(r != row) {
					this.globalLoad.add(r,-this.globalStiff.get(r, col)*value);
					this.globalStiff.set(r, col, 0.0);
				}
			}
		} else { //value=0
			for(int r=1;r<=nRow;r++) {
				if(r != row) {
					this.globalStiff.set(r, col, 0.0);
				}
			}
		}
		for(int c=1;c<=nCol;c++) {
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
	public void imposeDirichletCondition(VecMathFunc diri) {
		ElementList eList = mesh.getElementList();
		
		int nRow = this.globalStiff.getRowDim();
		int nCol = this.globalStiff.getColDim();
		FullMatrix fStiff = new FullMatrix(nRow,nCol);
		double [][]fsData = fStiff.getData();
		//System.out.println(this.globalStiff.getNonZeroNumber());
		Map<Integer,Map<Integer,Double>> sData = this.globalStiff.getAll();
		//int size = 0;
		for(Entry<Integer,Map<Integer,Double>> e1 : sData.entrySet()) {
			int r = e1.getKey();
			Map<Integer,Double> row = e1.getValue();
			for(Entry<Integer,Double> e2 : row.entrySet()) {
				int c = e2.getKey();
				Double v = e2.getValue();
				fsData[r-1][c-1] = v;
				//size++;
			}
		}
		//System.out.println(this.globalStiff.getNonZeroNumber()+"=="+size);
		
		Set<Integer> nodeDOFSet = new HashSet<Integer>();
		for(int ie=1;ie<=eList.size();ie++) {
			Element e = eList.at(ie);
			DOFList DOFs = e.getAllDOFList(DOFOrder.NEFV);
			for(int j=1;j<=DOFs.size();j++) {
				DOF dof = DOFs.at(j);
				GeoEntity ge = dof.getOwner();
				int nVVFCmpt = dof.getVVFComponent();
				MathFunc fdiri = diri.get(nVVFCmpt);
				if(ge instanceof Node) {
					//if(!nodeDOFSet.contains(dof.getGlobalIndex())) {
						Node n = (Node)ge;
						if(n.getNodeType(nVVFCmpt) == NodeType.Dirichlet) {
							Variable v = Variable.createFrom(fdiri, n, 0);
							setDirichlet(fStiff,dof.getGlobalIndex(),fdiri.apply(v));
						}
					//	nodeDOFSet.add(dof.getGlobalIndex());
					//}
				} else if(ge instanceof EdgeLocal) {
					//2D单元（面）其中的局部边上的自由度
					EdgeLocal edge = (EdgeLocal)ge;
					if(edge.getBorderType(nVVFCmpt) == NodeType.Dirichlet) {
						//TODO 以边的那个顶点取值？中点？
						//Variable v = Variable.createFrom(fdiri, ?, 0);
					}
					
				} else if(ge instanceof FaceLocal) {
					//3D单元（体）其中的局部面上的自由度
					FaceLocal face = (FaceLocal)ge;
					if(face.getBorderType(nVVFCmpt) == NodeType.Dirichlet) {
						//TODO
					}
				} else if(ge instanceof Edge) {
					//1D单元（线段）上的自由度，其Dirichlet边界用结点来计算推出，而不需要专门标记单元
					VertexList vs = ((GeoEntity2D) ge).getVertices();
					for(int k=1;k<=vs.size();k++) {
						Node n = vs.at(k).globalNode();
						if(NodeType.Dirichlet == n.getNodeType(nVVFCmpt)) {
							Variable v = Variable.createFrom(fdiri, n, 0);
							setDirichlet(fStiff,dof.getGlobalIndex(),fdiri.apply(v));
						}
					}
				} else if(ge instanceof Face) {
					//2D单元（面）上的自由度，其Dirichlet边界用结点来计算推出，而不需要专门标记单元
					
					VertexList vs = ((GeoEntity2D) ge).getVertices();
					for(int k=1;k<=vs.size();k++) {
						Node n = vs.at(k).globalNode();
						if(NodeType.Dirichlet == n.getNodeType(nVVFCmpt)) {
							Variable v = Variable.createFrom(fdiri, n, 0);
							setDirichlet(fStiff,dof.getGlobalIndex(),fdiri.apply(v));
						}
					}
				} else if(ge instanceof Volume) {
					//3D单元（体）上的自由度，其Dirichlet边界用结点来计算推出，而不需要专门标记单元
					VertexList vs = ((GeoEntity3D) ge).getVertices();
					for(int k=1;k<=vs.size();k++) {
						Node n = vs.at(k).globalNode();
						if(NodeType.Dirichlet == n.getNodeType(nVVFCmpt)) {
							Variable v = Variable.createFrom(fdiri, n, 0);
							setDirichlet(fStiff,dof.getGlobalIndex(),fdiri.apply(v));
						}
					}
				}
			}
		}
		this.globalStiff.clearData();
		long begin = System.currentTimeMillis();
		for(int i=nRow; --i>=0; ) {
		//for(int i=0; i<nRow; i++) {
			double[] _fsDatai = fsData[i];
			for(int j=nCol; --j>=0; ) {
			//for(int j=0; j<nCol; j++) {
				double v = _fsDatai[j];
				if(Math.abs(v) > Matrix.zeroEps) {
					this.globalStiff.set(i+1, j+1, v);
				}
			}
		}
		long end = System.currentTimeMillis();
		System.out.println("time="+(end-begin)+"ms");
	}		
}
