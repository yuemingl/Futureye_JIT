package edu.uta.futureye.core;

import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import edu.uta.futureye.algebra.Matrix;
import edu.uta.futureye.algebra.Vector;
import edu.uta.futureye.core.intf.WeakForm;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.util.DOFList;
import edu.uta.futureye.util.ElementList;
import edu.uta.futureye.util.NodeList;
import edu.uta.futureye.util.PairElementMatrix;

public class Assembler {
	protected Mesh mesh;
	protected WeakForm weakForm;
	protected Matrix globalMatrix;
	protected Vector loadVector;

	public Assembler(Mesh mesh, WeakForm weakForm) {
		this.mesh = mesh;
		this.weakForm = weakForm;
		
		int dim = mesh.nodeList.size();
		globalMatrix = new Matrix(dim,dim);
		loadVector = new Vector(dim);

	}
	
	public Matrix getStiffnessMatrix() {
		
		ElementList eleList = mesh.getElementList();
		for(int i=1;i<=eleList.size();i++) {
			Element e = eleList.at(i);
			e.adjustVerticeToCounterClockwise();
			getLocalMatrix(e);
		}
		procHangingNode(mesh);
		
		return globalMatrix;
	}
	
	//二维：刚度矩阵增加hanging node约束系数
	// nh - 0.5*n1 - 0.5*n2 = 0
	public void procHangingNode(Mesh mesh) {
		
		for(int i=1;i<=mesh.getNodeList().size();i++) {
			Node node = mesh.getNodeList().at(i);
			if(node instanceof NodeRefined) {
				NodeRefined nRefined = (NodeRefined)node;
				if(nRefined.isHangingNode()) {
					globalMatrix.set(nRefined.globalIndex, nRefined.globalIndex, 1.0);
					globalMatrix.set(nRefined.globalIndex,
							nRefined.constrainNodes.at(1).globalIndex,-0.5);
					globalMatrix.set(nRefined.globalIndex,
							nRefined.constrainNodes.at(2).globalIndex,-0.5);
				}
			}
		}
	}
	
	public Vector getLoadVector() {
		ElementList eleList = mesh.getElementList();
		for(int i=1;i<=eleList.size();i++) {
			Element e = eleList.at(i);
			plusToLoadVector(e);
		}
		return loadVector;
	}
	
	public void plusToLoadVector(Element e) {
		int nNode = e.nodes.size();
		
		e.updateJacobinLinear2D();
		for(int j=1;j<=nNode;j++) {
			DOFList dofListJ = e.getDOFList(j);
			int nDOF = dofListJ.size();
			for(int k=1;k<=nDOF;k++) {
				DOF dofJ = dofListJ.at(k);
				weakForm.setShapeFunction(
							null,0,
							dofJ.getShapeFunction(),dofJ.localIndex
						);
				double val = weakForm.rightHandSide(e, WeakForm.ItemType.Domain).value(null);
				//System.out.println(val);
				loadVector.plusValue(dofJ.globalIndex, val);
			}
		}
		
		if(e.isBorderElement()) {
			ElementList beList = e.getSubElements();
			//System.out.println(">>Border:"+beList);
			for(int ibeList=1;ibeList<=beList.size();ibeList++) {
				Element be = beList.at(ibeList);
				Edge edge = (Edge)be.stemGeoEntity;
				if(edge.getBorderNodeType() == NodeType.Neumann ||
						edge.getBorderNodeType() == NodeType.Robin) {
					be.updateJacobinLinear1D();
					nNode = be.nodes.size();
					for(int j=1;j<=nNode;j++) {
						DOFList dofListJ = be.getDOFList(j);
						int nDOF = dofListJ.size();
						for(int k=1;k<=nDOF;k++) {
							DOF dofJ = dofListJ.at(k);
							//TODO
							weakForm.setShapeFunction(
									null,0,
									dofJ.getShapeFunction(),dofJ.localIndex
									);
							
							//lym leftHandSide -> rightHandSide
							Function fun = weakForm.rightHandSide(be, WeakForm.ItemType.Border);
							if(fun != null) { //TODO 如果返回null，说明此项调用对于该弱形式不适用，不进行载荷向量的累加
							loadVector.plusValue(dofJ.globalIndex, 
									fun.value(null));
							//loadVector.set(dofJ.globalIndex, 0);
							}
						}
					}
				}
			}
		}
	}
	
	
	public Matrix getLocalMatrix(Element e) {
		//System.out.println("Assembling:"+e);
		Matrix localMatrix = null;
		List<PairElementMatrix> pemList = weakForm.associateElement(e);
		
		int dofDim = 0;
		if(pemList == null) {
			dofDim = e.getTotalNumberOfDOF();
			localMatrix = new Matrix(dofDim,dofDim);
		} else {
			for(PairElementMatrix pem : pemList) {
				plusToGlobalMatrix(pem.element, pem.localMatrix);
			}
			return null; //TODO ???
		}
		
		e.updateJacobinLinear2D();
		int nNode = e.nodes.size();
		//结点双循环
		for(int i=1;i<=nNode;i++) {
			for(int j=1;j<=nNode;j++) {
				DOFList dofListI = e.getDOFList(i);
				DOFList dofListJ = e.getDOFList(j);
				
				//结点上自由度双循环
				//不同结点上的自由度个数不一定相同(e.g. adaptive finite element)
				int nDOF1 = dofListI.size();
				int nDOF2 = dofListJ.size();
				for(int k1=1;k1<=nDOF1;k1++) {
					DOF dofI = dofListI.at(k1);
					for(int k2=1;k2<=nDOF2;k2++) {
						DOF dofJ = dofListJ.at(k2);
						//TODO
						int lRow = dofI.localIndex;
						int lCol = dofJ.localIndex;
						weakForm.setShapeFunction(
								dofI.getShapeFunction(),dofI.localIndex,
								dofJ.getShapeFunction(),dofJ.localIndex
								);
						double val = weakForm.leftHandSide(e, WeakForm.ItemType.Domain).value(null);
						//System.out.println(e+" "+val);
						localMatrix.set(lRow, lCol, val);
					}
				}
			}
		}
		plusToGlobalMatrix(e, localMatrix);
		
		if(e.isBorderElement()) {
			ElementList beList = e.getSubElements();
			//System.out.println(">>Border:"+beList);	
			for(int ibeList=1;ibeList<=beList.size();ibeList++) {
			
				Element be = beList.at(ibeList);
				dofDim = be.getTotalNumberOfDOF();
				Edge edge = (Edge)be.stemGeoEntity;
				if(edge.getBorderNodeType() == NodeType.Robin) {
					be.updateJacobinLinear1D();
					nNode = be.nodes.size();
					localMatrix = new Matrix(dofDim,dofDim);
					for(int i=1;i<=nNode;i++) {
						for(int j=1;j<=nNode;j++) {
							DOFList dofListI = be.getDOFList(i);
							DOFList dofListJ = be.getDOFList(j);
							
							int nDOF1 = dofListI.size();
							int nDOF2 = dofListJ.size();
							for(int k1=1;k1<=nDOF1;k1++) {
								DOF dofI = dofListI.at(k1);
								for(int k2=1;k2<=nDOF2;k2++) {
									DOF dofJ = dofListJ.at(k2);
									//TODO
									int lRow = dofI.localIndex;
									int lCol = dofJ.localIndex;
									weakForm.setShapeFunction(
											dofI.getShapeFunction(),dofI.localIndex,
											dofJ.getShapeFunction(),dofJ.localIndex
											);
									Function fun = weakForm.leftHandSide(be, WeakForm.ItemType.Border);
									if(fun != null) { //TODO 如果返回null，说明此项调用对于该弱形式不适用，不进行刚度矩阵的累加
										double rlt = fun.value(null);
										//System.out.println(rlt);
										localMatrix.plusValue(lRow, lCol, rlt);
									}
								}
							}
						}
					}
					plusToGlobalMatrix(be, localMatrix);
				}
			}
		}
		return localMatrix;
	}
	
	
	public void plusToGlobalMatrix(Element e, Matrix local) {
		
		//local.print();
		
		Map<Integer,Map<Integer,Double>> m = local.getAll();
		for(Entry<Integer,Map<Integer,Double>> er : m.entrySet()) {
			int nRow = er.getKey();
			int ngRow = e.local2GlobalDOFIndex(nRow);
			
			//System.out.println(nRow+"->"+ngRow);
			
			Map<Integer,Double> row = er.getValue();
			for(Entry<Integer,Double> ec : row.entrySet()) {
				int nCol = ec.getKey();
				double val = ec.getValue();
				int ngCol = e.local2GlobalDOFIndex(nCol);
				
				this.globalMatrix.plusValue(ngRow, ngCol, val);
			}
		}
	}
	
	public void imposeDirichletCondition(Function diri) {
		NodeList nList = mesh.nodeList;
		for(int i=1;i<=nList.size();i++) {
			Node n = nList.at(i);
			if(n.getNodeType() == NodeType.Dirichlet) {
				Variable v = Variable.createFrom(diri, n, n.globalIndex);
				this.globalMatrix.set(n.globalIndex, n.globalIndex, 1.0);
				this.loadVector.set(n.globalIndex, diri.value(v));
				for(int j=1;j<=this.globalMatrix.getRowDim();j++) {
					if(j!=n.globalIndex) {
						//this.globalMatrix.set(j, n.globalIndex, 0.0);
						this.globalMatrix.set(n.globalIndex, j, 0.0);
					}
				}
			}
		}
	}
	
}
