package edu.uta.futureye.lib.assembler;

import edu.uta.futureye.algebra.SparseMatrixRowMajor;
import edu.uta.futureye.algebra.SparseVectorHashMap;
import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.core.DOF;
import edu.uta.futureye.core.DOFOrder;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.function.operator.FOIntegrate;
import edu.uta.futureye.lib.weakform.WeakForm;
import edu.uta.futureye.util.container.DOFList;
import edu.uta.futureye.util.container.ElementList;

public class Assembler {
	double[][] A; // local stiff matrix
	double[] b;   // local load vector
	double[] params;
	WeakForm domainWF;
	WeakForm boundaryWF;
	int nDOFs;
	
	Matrix gA; // global stiff matrix
	Vector gb; // global load vector

	public Assembler(WeakForm domainWeakForm) {
		this.domainWF = domainWeakForm;
		nDOFs = domainWF.getFiniteElement().getNumberOfDOFs();
		A = new double[nDOFs][nDOFs];
		b = new double[nDOFs];
		params = new double[domainWF.getFiniteElement().getArgsOrder().length];
	}
	
	/**
	 * 
	 * @param domainWeakForm
	 * @param boundaryWeakForm
	 */
	public Assembler(WeakForm domainWeakForm, WeakForm boundaryWeakForm) {
		this.domainWF = domainWeakForm;
		this.boundaryWF = boundaryWeakForm;
		nDOFs = domainWF.getFiniteElement().getNumberOfDOFs();
		A = new double[nDOFs][nDOFs];
		b = new double[nDOFs];
		params = new double[domainWF.getFiniteElement().getArgsOrder().length];
	}
	
	/**
	 * Assemble on a give element
	 * @param e
	 */
	public void assembleLocal(Element e) {
		//e.adjustVerticeToCounterClockwise();

		double[] coords = e.getNodeCoords();
		System.arraycopy(coords, 0, params, 0, coords.length);

		domainWF.getCompiledJac().apply(params);

		if(domainWF.getFiniteElement().getNumberOfDOFs() == 3) {
			for(int j=0; j<nDOFs; j++) {
				for(int i=0; i<nDOFs; i++) {
					A[j][i] = FOIntegrate.intOnTriangleRefElement(domainWF.getCompiledLHS()[j][i], 
							params, coords.length, 2);//2=80.839 3=80.966, 4=80.967
				}
				b[j] = FOIntegrate.intOnTriangleRefElement(domainWF.getCompiledRHS()[j], 
						params, coords.length, 2);
			}
		} else if(domainWF.getFiniteElement().getNumberOfDOFs() == 4) {
			for(int j=0; j<nDOFs; j++) {
				for(int i=0; i<nDOFs; i++) {
					A[j][i] = FOIntegrate.intOnRectangleRefElement(domainWF.getCompiledLHS()[j][i], 
							params, coords.length, 2);
				}
				b[j] = FOIntegrate.intOnRectangleRefElement(domainWF.getCompiledRHS()[j], 
						params, coords.length, 2);
			}
		}
	}
	
	/**
	 * Assemble stiff matrix and load vector on a given mesh
	 * @param mesh
	 */
	public void assembleGlobal(Mesh mesh) {
		int dim = mesh.getNodeList().size();
		gA = new SparseMatrixRowMajor(dim,dim);
		gb = new SparseVectorHashMap(dim);
		assembleGlobal(mesh, gA, gb);
	}
	
	/**
	 * Assemble stiff matrix and load vector on a given mesh
	 * into parameter stiff and load.
	 * 
	 * Several assemblers can be chained by using this method
	 * to assemble stiff matrix and load vector
	 * 
	 * @param mesh
	 * @param stiff
	 * @param load
	 */
	public void assembleGlobal(Mesh mesh, Matrix stiff, Vector load) {
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
					stiff.add(nGlobalRow, nGlobalCol, A[j][i]);
				}
				//Local load vector
				load.add(nGlobalRow, b[j]);
			}
		}
		//update gA and gb
		this.gA = stiff;
		this.gb = load;
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
