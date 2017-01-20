package edu.uta.futureye.lib.element;

import edu.uta.futureye.algebra.SparseMatrixRowMajor;
import edu.uta.futureye.algebra.SparseVectorHashMap;
import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.SparseMatrix;
import edu.uta.futureye.algebra.intf.SparseVector;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.core.DOF;
import edu.uta.futureye.core.DOFOrder;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.function.operator.FOIntegrate;
import edu.uta.futureye.util.container.DOFList;
import edu.uta.futureye.util.container.ElementList;

public class AssemblerJIT {
	double[][] A; // local stiff matrix
	double[] b;   // local load vector
	double[] params;
	WeakFormJIT wf;
	int nDOFs;
	
	SparseMatrix gA; // global stiff matrix
	SparseVector gb; // global load vector

	public AssemblerJIT(WeakFormJIT wf) {
		this.wf = wf;
		nDOFs = wf.getFiniteElement().getNumberOfDOFs();
		A = new double[nDOFs][nDOFs];
		b = new double[nDOFs];
		params = new double[wf.getFiniteElement().getArgsOrder().length];
	}
	
	/**
	 * Assemble on a give element
	 * @param e
	 */
	public void assembleLocal(Element e) {
		//e.adjustVerticeToCounterClockwise();

		double[] coords = e.getNodeCoords();
		System.arraycopy(coords, 0, params, 0, coords.length);

		wf.getCompiledJac().apply(params);

		for(int j=0; j<nDOFs; j++) {
			for(int i=0; i<nDOFs; i++) {
				A[j][i] = FOIntegrate.intOnTriangleRefElement(wf.getCompiledLHS()[j][i], 
						params, coords.length, 2);//2=80.839 3=80.966, 4=80.967
			}
			b[j] = FOIntegrate.intOnTriangleRefElement(wf.getCompiledRHS()[j], 
					params, coords.length, 2);
		}
	}
	
	/**
	 * Assemble on a given mesh
	 * @param mesh
	 */
	public void assembleGlobal(Mesh mesh) {
		int dim = mesh.getNodeList().size();
		gA = new SparseMatrixRowMajor(dim,dim);
		gb = new SparseVectorHashMap(dim);
		ElementList eList = mesh.getElementList();
		
		for(Element e : eList) {
			assembleLocal(e);
			DOFList DOFs = e.getAllDOFList(DOFOrder.NEFV);
			for(int j=0;j<nDOFs;j++) {
				DOF dofI = DOFs.at(j+1);
				int nGlobalRow = dofI.getGlobalIndex();
				for(int i=0;i<nDOFs;i++) {
					DOF dofJ = DOFs.at(i+1);
					int nGlobalCol = dofJ.getGlobalIndex();
					gA.add(nGlobalRow, nGlobalCol, A[j][i]);
				}
				//Local load vector
				gb.add(nGlobalRow, b[j]);
			}
		}
	}
	
	public double[][] getLocalStiffMatrix() {
		return A;
	}
	
	public double[] getLocalLoadVector() {
		return b;
	}
	
	public Matrix getGlobalStiffMatrix() {
		return gA;
	}

	public Vector getGlobalLoadVector() {
		return gb;
	}
}
