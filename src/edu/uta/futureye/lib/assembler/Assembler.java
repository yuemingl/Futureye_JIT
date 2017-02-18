/**
 * Copyright (c) 2010, nkliuyueming@gmail.com. All rights reserved.
 * 
 * 
 */
package edu.uta.futureye.lib.assembler;

import edu.uta.futureye.algebra.SparseMatrixRowMajor;
import edu.uta.futureye.algebra.SparseVectorHashMap;
import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.intf.FiniteElement;
import edu.uta.futureye.lib.weakform.WeakForm;

public class Assembler {
	Mesh mesh;
	BasicAssembler assembler;
	
	Matrix gA;
	Vector gb;
	
	Assembler parentAssembler;
	
	public Assembler(Mesh mesh, WeakForm weakForm) {
		this.mesh = mesh;
		this.assembler = new BasicAssembler(mesh, weakForm);
	}

	/**
	 * Several assemblers can be chained by using this method
	 * to assemble stiff matrix and load vector
	 * @param weakForm
	 */
	public Assembler(Assembler parent, WeakForm weakForm) {
		this.parentAssembler = parent;
		this.assembler = new BasicAssembler(mesh, weakForm);
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
	public void assembleGlobal() {
		int dim = this.assembler.weakForm.getFiniteElement().getTotalNumberOfDOFs(mesh);
		if(null == this.parentAssembler) {
			gA = new SparseMatrixRowMajor(dim,dim);
			gb = new SparseVectorHashMap(dim);
			assembleGlobal(gA, gb);
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
	public void assembleGlobal(Matrix stiff, Vector load) {
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
		FiniteElement fe = this.assembler.weakForm.getFiniteElement();
		int nDOFs = fe.getNumberOfDOFs();

		// Get local-global indexing
		for(int j=0;j<nDOFs;j++) {
			int globalIdxJ = fe.getGlobalIndex(mesh, e, j+1);
			for(int i=0;i<nDOFs;i++) {
				int globalIdxI = fe.getGlobalIndex(mesh, e, i+1);
				stiff.add(globalIdxJ, globalIdxI, this.assembler.A[j][i]);
			}
			load.add(globalIdxJ, this.assembler.b[j]);
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
