package edu.uta.futureye.lib.assembler;

import edu.uta.futureye.algebra.SparseBlockMatrix;
import edu.uta.futureye.algebra.SparseBlockVector;
import edu.uta.futureye.algebra.SparseMatrix;
import edu.uta.futureye.algebra.SparseVector;
import edu.uta.futureye.algebra.intf.BlockMatrix;
import edu.uta.futureye.algebra.intf.BlockVector;
import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.intf.Assembler;
import edu.uta.futureye.core.intf.WeakForm;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.function.intf.VectorFunction;
import edu.uta.futureye.util.container.ElementList;

public class AssemblerMixedLaplace implements Assembler {
	protected Mesh mesh;
	protected WeakForm weakForm;
	protected BlockMatrix globalStiff;
	protected BlockVector globalLoad;

	public AssemblerMixedLaplace(Mesh mesh, WeakForm weakForm) {
		this.mesh = mesh;
		this.weakForm = weakForm;
		
		int edgeDOF = mesh.getEdgeList().size();
		int elementDOF = mesh.getEdgeList().size();
		
		globalStiff = new SparseBlockMatrix(2,2);
		globalStiff.setBlock(1, 1, new SparseMatrix(edgeDOF,edgeDOF));
		globalStiff.setBlock(1, 2, new SparseMatrix(edgeDOF,elementDOF));
		globalStiff.setBlock(2, 1, new SparseMatrix(elementDOF,edgeDOF));
		globalStiff.setBlock(2, 2, new SparseMatrix(elementDOF,elementDOF));
		
		globalLoad = new SparseBlockVector(2);
		globalLoad.setBlock(1, new SparseVector(edgeDOF));
		globalLoad.setBlock(2, new SparseVector(elementDOF));
		
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
	public Vector getLoadVector() {
		return globalLoad;
	}

	@Override
	public Matrix getStiffnessMatrix() {
		return globalStiff;
	}

	@Override
	public void imposeDirichletCondition(Function diri) {
		//不需要任何处理
		throw new UnsupportedOperationException();
	}

	@Override
	public void imposeDirichletCondition(VectorFunction diri) {
		throw new UnsupportedOperationException();
	}
}
