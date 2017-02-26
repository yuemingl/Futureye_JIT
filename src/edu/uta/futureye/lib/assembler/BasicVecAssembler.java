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
import edu.uta.futureye.core.intf.VecFiniteElement;
import edu.uta.futureye.function.operator.FOIntegrate;
import edu.uta.futureye.lib.weakform.VecWeakForm;
import edu.uta.futureye.util.container.ElementList;

public class BasicVecAssembler {
	public Mesh mesh;
	public VecWeakForm weakForm;
	public double[][] A; // domain local stiff matrix
	public double[] b;   // domain local load vector
	double[] params;
	int nDOFs;
	
	Matrix gA; // global stiff matrix
	Vector gb; // global load vector

	public BasicVecAssembler(Mesh mesh, VecWeakForm weakForm) {
		this.mesh = mesh;
		this.weakForm = weakForm;
		this.nDOFs = weakForm.getFiniteElement().getNumberOfDOFs();
		this.A = new double[nDOFs][nDOFs];
		this.b = new double[nDOFs];
		this.params = new double[weakForm.getFiniteElement().getArgsOrder().length];
	}
	
	/**
	 * Assemble local stiff matrix and load vector on a give element
	 * @param e
	 */
	public void assembleLocal(Element e) {
		e.adjustVerticeToCounterClockwise();
		
		double[] coords = e.getNodeCoords();
		System.arraycopy(coords, 0, params, 0, coords.length);

		weakForm.getCompiledJac().apply(params);

		if(e.vertices().size() == 2) {
			for(int j=0;j<nDOFs;j++) {
				for(int i=0;i<nDOFs;i++) {
					A[j][i] = FOIntegrate.intOnLinearRefElement(weakForm.getCompiledLHS()[j][i], 
							new AssembleParam(e, i+1, j+1), params, coords.length, 5);
				}
				b[j] = FOIntegrate.intOnLinearRefElement(weakForm.getCompiledRHS()[j], 
						new AssembleParam(e, -1, j+1), params, coords.length, 5);
			}
		} else if(e.vertices().size() == 3) {
			
			for(int j=0; j<nDOFs; j++) {
				for(int i=0; i<nDOFs; i++) {
					
					//this.weakForm.getFiniteElement().getGeoEntity(e, i);
					
					A[j][i] = FOIntegrate.intOnTriangleRefElement(weakForm.getCompiledLHS()[j][i], 
							new AssembleParam(e, i+1, j+1), params, coords.length, 5); //TODO use high order for high order element???
					//Laplace Test: 2=80.839 3=80.966, 4=80.967
				}
				b[j] = FOIntegrate.intOnTriangleRefElement(weakForm.getCompiledRHS()[j], 
						new AssembleParam(e, -1, j+1), params, coords.length, 5);
			}
		} else if(e.vertices().size() == 4) {
			for(int j=0; j<nDOFs; j++) {
				for(int i=0; i<nDOFs; i++) {
					A[j][i] = FOIntegrate.intOnRectangleRefElement(weakForm.getCompiledLHS()[j][i], 
							new AssembleParam(e, i+1, j+1), params, coords.length, 5);
				}
				b[j] = FOIntegrate.intOnRectangleRefElement(weakForm.getCompiledRHS()[j], 
						new AssembleParam(e, -1, j+1), params, coords.length, 5);
			}
		}
		
	}
	
	/**
	 * Assemble global stiff matrix and load vector on a given mesh
	 * new matrix and vector are allocated. Use <tt>getGlobalStiffMatrix()</tt> and
	 * <tt>getGlobalLoadVector()</tt> to access them.
	 */
	public void assembleGlobal() {
		int dim = this.weakForm.getFiniteElement().getTotalNumberOfDOFs(mesh);
		
		gA = new SparseMatrixRowMajor(dim,dim);
		gb = new SparseVectorHashMap(dim);
		
		assembleGlobal(gA, gb);
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
		ElementList eList = mesh.getElementList();
		VecFiniteElement fe = this.weakForm.getFiniteElement();
		for(Element e : eList) {
			assembleLocal(e);

			for(int j=0;j<nDOFs;j++) {
				int nGlobalRow = fe.getGlobalIndex(mesh, e, j+1);
				for(int i=0;i<nDOFs;i++) {
					int nGlobalCol = fe.getGlobalIndex(mesh, e, i+1);
//					if(nGlobalRow==917) {
//						System.out.println("(i,j)=("+i+","+j+"); global=("+nGlobalRow+","+nGlobalCol+"); value="+A[j][i]);
//					}
					//if(Math.abs(A[j][i])>1e-3) {
						stiff.add(nGlobalRow, nGlobalCol, A[j][i]);
					//}
				}
				//Local load vector
				load.add(nGlobalRow, b[j]);
			}
		}
		//update gA and gb
		this.gA = stiff;
		this.gb = load;
	}

	public void assembleGlobal(Element e, Matrix stiff, Vector load) {
		// Assemble locally
		assembleLocal(e);
		VecFiniteElement fe = this.weakForm.getFiniteElement();
		int nDOFs = fe.getNumberOfDOFs();

		// Get local-global indexing
		for(int j=0;j<nDOFs;j++) {
			int globalIdxJ = fe.getGlobalIndex(mesh, e, j+1);
			for(int i=0;i<nDOFs;i++) {
				int globalIdxI = fe.getGlobalIndex(mesh, e, i+1);
				stiff.add(globalIdxJ, globalIdxI, A[j][i]);
			}
			load.add(globalIdxJ, b[j]);
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
