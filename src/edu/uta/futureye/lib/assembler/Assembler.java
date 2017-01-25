package edu.uta.futureye.lib.assembler;

import edu.uta.futureye.algebra.SparseMatrixRowMajor;
import edu.uta.futureye.algebra.SparseVectorHashMap;
import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.core.DOF;
import edu.uta.futureye.core.DOFOrder;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.lib.weakform.WeakForm;
import edu.uta.futureye.util.container.DOFList;

public class Assembler {
	BasicAssembler assembler;
	
	Matrix gA;
	Vector gb;
	
	Assembler parentAssembler;
	
	public Assembler(WeakForm weakForm) {
		this.assembler = new BasicAssembler(weakForm);
	}

	/**
	 * Several assemblers can be chained by using this method
	 * to assemble stiff matrix and load vector
	 * @param weakForm
	 */
	public Assembler(Assembler parent, WeakForm weakForm) {
		this.parentAssembler = parent;
		this.assembler = new BasicAssembler(weakForm);
	}
	
	/**
	 * Assemble local stiff and load on a give element
	 * @param e
	 */
	public void assembleLocal(Element e) {
		// Assemble on domain element
		assembler.assembleLocal(e);
	}
	
	
	/**
	 * Assemble stiff matrix and load vector on a given mesh
	 * @param mesh
	 */
	public void assembleGlobal(Element e) {
		if(null != this.parentAssembler) {
			assembleGlobal(e, 
					this.parentAssembler.getGlobalStiffMatrix(), 
					this.parentAssembler.getGlobalLoadVector());
		} else {
			throw new RuntimeException("Call assembleGlobal(Mesh mesh) first!");
		}
	}

	/**
	 * Assemble stiff matrix and load vector on a given mesh
	 * @param mesh
	 */
	public void assembleGlobal(Mesh mesh) {
		int dim = mesh.getNodeList().size();
		if(null == this.parentAssembler) {
			gA = new SparseMatrixRowMajor(dim,dim);
			gb = new SparseVectorHashMap(dim);
			assembleGlobal(mesh, gA, gb);
		} else {
			throw new RuntimeException("Call assembleGlobal(Mesh mesh) in root assembler only!");
		}
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
		for(Element e : mesh.getElementList()) {
			this.assembleGlobal(e, stiff, load);
		}
	}
	
	/**
	 * 
	 * @param e
	 * @param stiff
	 * @param load
	 */
	public void assembleGlobal(Element e, Matrix stiff, Vector load) {
		// Assemble locally
		assembleLocal(e);

		// Get local-global indexing
		DOFList DOFs = e.getAllDOFList(DOFOrder.NEFV);
		for(int j=0;j<DOFs.size();j++) {
			DOF dofJ = DOFs.at(j+1);
			for(int i=0;i<DOFs.size();i++) {
				DOF dofI = DOFs.at(i+1);
				stiff.add(dofJ.getGlobalIndex(), dofI.getGlobalIndex(), this.assembler.A[j][i]);
			}
			load.add(dofJ.getGlobalIndex(), this.assembler.b[j]);
		}
		//update gA and gb
		this.gA = stiff;
		this.gb = load;
	}

	public double[][] getLocalStiffMatrix() {
		return this.assembler.getLocalStiffMatrix();
	}
	
	public double[] getLocalLoadVector() {
		return this.assembler.getLocalLoadVector();
	}

	public Matrix getGlobalStiffMatrix() {
		return gA;
	}

	public Vector getGlobalLoadVector() {
		return gb;
	}
}
