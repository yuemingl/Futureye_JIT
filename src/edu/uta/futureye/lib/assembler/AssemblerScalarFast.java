package edu.uta.futureye.lib.assembler;

import edu.uta.futureye.algebra.SparseMatrixRowMajor;
import edu.uta.futureye.algebra.SparseVectorHashMap;
import edu.uta.futureye.algebra.intf.SparseMatrix;
import edu.uta.futureye.algebra.intf.SparseVector;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.core.NodeRefined;
import edu.uta.futureye.core.NodeType;
import edu.uta.futureye.core.intf.AssemblerOld;
import edu.uta.futureye.core.intf.WeakFormOld;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.function.intf.VecMathFunc;
import edu.uta.futureye.util.container.ElementList;
import edu.uta.futureye.util.container.NodeList;

public class AssemblerScalarFast implements AssemblerOld{
	protected Mesh mesh;
	protected WeakFormOld weakForm;
	protected SparseMatrix globalStiff;
	protected SparseVector globalLoad;

	public AssemblerScalarFast(Mesh mesh, WeakFormOld weakForm) {
		this.mesh = mesh;
		this.weakForm = weakForm;
		
		int dim = mesh.getNodeList().size();
		globalStiff = new SparseMatrixRowMajor(dim,dim);
		globalLoad = new SparseVectorHashMap(dim);

	}
	
	@Override
	public void assemble() {
		ElementList eList = mesh.getElementList();
		int nEle = eList.size();
		for(int i=1; i<=nEle; i++) {
			eList.at(i).adjustVerticeToCounterClockwise();
			weakForm.assembleElement(eList.at(i), 
					globalStiff, globalLoad);
			if(i%3000==0)
				System.out.println("Assemble..."+
						String.format("%.0f%%", 100.0*i/nEle));
		}
		procHangingNode(mesh);
		return;
	}
	
	@Override
	public SparseVector getLoadVector() {
		return globalLoad;
	}

	@Override
	public SparseMatrix getStiffnessMatrix() {
		return globalStiff;
	}
	
	@Override
	public void imposeDirichletCondition(MathFunc diri) {
		NodeList nList = mesh.getNodeList();
		for(int i=1;i<=nList.size();i++) {
			Node n = nList.at(i);
			if(n.getNodeType() == NodeType.Dirichlet) {
				Variable v = Variable.createFrom(diri, n, n.globalIndex);
				this.globalStiff.set(n.globalIndex, n.globalIndex, 1.0);
				double val = diri.apply(v);
				this.globalLoad.set(n.globalIndex, val);
				for(int j=1;j<=this.globalStiff.getRowDim();j++) {
					if(j!=n.globalIndex) {
						//TODO 行列都需要置零
						this.globalLoad.add(j, -this.globalStiff.get(j, n.globalIndex)*val);
						this.globalStiff.set(j, n.globalIndex, 0.0);
						this.globalStiff.set(n.globalIndex, j, 0.0);
					}
				}
			}
		}
	}
	
	//二维：刚度矩阵增加hanging node约束系数
	// nh - 0.5*n1 - 0.5*n2 = 0
	public void procHangingNode(Mesh mesh) {
		
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
	
	@Override
	public void imposeDirichletCondition(VecMathFunc diri) {
		throw new UnsupportedOperationException();
	}	
}
