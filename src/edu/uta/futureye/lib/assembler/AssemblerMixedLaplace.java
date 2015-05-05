package edu.uta.futureye.lib.assembler;

import edu.uta.futureye.algebra.SparseBlockMatrix;
import edu.uta.futureye.algebra.SparseBlockVector;
import edu.uta.futureye.algebra.SparseMatrixRowMajor;
import edu.uta.futureye.algebra.SparseVectorHashMap;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.intf.Assembler;
import edu.uta.futureye.core.intf.WeakForm;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.function.intf.VectorFunction;
import edu.uta.futureye.util.container.ElementList;

public class AssemblerMixedLaplace implements Assembler {
	protected Mesh mesh;
	protected WeakForm weakForm;
	protected SparseBlockMatrix globalStiff;
	protected SparseBlockVector globalLoad;

	public AssemblerMixedLaplace(Mesh mesh, WeakForm weakForm) {
		this.mesh = mesh;
		this.weakForm = weakForm;
		
		int edgeDOF = mesh.getEdgeList().size();
		int elementDOF = mesh.getEdgeList().size();
		
		globalStiff = new SparseBlockMatrix(2,2);
		globalStiff.setBlock(1, 1, 
				new SparseMatrixRowMajor(edgeDOF,edgeDOF));
		globalStiff.setBlock(1, 2, 
				new SparseMatrixRowMajor(edgeDOF,elementDOF));
		globalStiff.setBlock(2, 1, 
				new SparseMatrixRowMajor(elementDOF,edgeDOF));
		globalStiff.setBlock(2, 2, 
				new SparseMatrixRowMajor(elementDOF,elementDOF));
		
		globalLoad = new SparseBlockVector(2);
		globalLoad.setBlock(1, new SparseVectorHashMap(edgeDOF));
		globalLoad.setBlock(2, new SparseVectorHashMap(elementDOF));
		
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
		return;
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
		//不需要任何处理
		throw new UnsupportedOperationException();
	}

	@Override
	public void imposeDirichletCondition(VectorFunction diri) {
		throw new UnsupportedOperationException();
	}
}
