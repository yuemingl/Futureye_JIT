package edu.uta.futureye.lib.assembler;

import edu.uta.futureye.algebra.SparseMatrixRowMajor;
import edu.uta.futureye.algebra.SparseVectorHashMap;
import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.bytecode.CompiledFunc;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.NodeType;
import edu.uta.futureye.core.intf.FiniteElement;
import edu.uta.futureye.function.operator.FOIntegrate;
import edu.uta.futureye.lib.weakform.WeakForm;
import edu.uta.futureye.util.container.ElementList;

public class DomainBoundaryAssemblerRaw {
	WeakForm domainWF;
	double[][] A; // domain local stiff matrix
	double[] b;   // domain local load vector
	double[] params;
	int nDOFs;
	
	WeakForm boundaryWF;
	double[][] beA; // boundary local stiff matrix
	double[] beb;   // boundary local load vector
	double[] beParams;
	int nBeDOFs;
	
	Matrix gA; // global stiff matrix
	Vector gb; // global load vector

	
	/**
	 * Construct an assembler for domain weak form and boundary weak form assemble
	 * 
	 * @param domainWeakForm
	 * @param boundaryWeakForm
	 */
	public DomainBoundaryAssemblerRaw(WeakForm domainWeakForm, WeakForm boundaryWeakForm) {
		this.domainWF = domainWeakForm;
		this.boundaryWF = boundaryWeakForm;
		
		nDOFs = domainWF.getFiniteElement().getNumberOfDOFs();
		A = new double[nDOFs][nDOFs];
		b = new double[nDOFs];
		params = new double[domainWF.getFiniteElement().getArgsOrder().length];
		
		if(null != boundaryWF) {
			nBeDOFs = this.boundaryWF.getFiniteElement().getNumberOfDOFs();
			beA = new double[nBeDOFs][nBeDOFs];
			beb = new double[nBeDOFs];
			beParams = new double[boundaryWF.getFiniteElement().getArgsOrder().length];
		}
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

		if(domainWF.getFiniteElement().getNumberOfDOFs() == 2) {
		for(int j=0;j<nDOFs;j++) {
			for(int i=0;i<nDOFs;i++) {
				A[j][i] = FOIntegrate.intOnLinearRefElement(domainWF.getCompiledLHS()[j][i], 
						new AssembleParam(e, i+1, j+1), params, coords.length, 5);
			}
			b[j] = FOIntegrate.intOnLinearRefElement(domainWF.getCompiledRHS()[j], 
					new AssembleParam(e, -1, j+1), params, coords.length, 5);
		}
		} else if(domainWF.getFiniteElement().getNumberOfDOFs() == 3) {
			for(int j=0; j<nDOFs; j++) {
				for(int i=0; i<nDOFs; i++) {
					A[j][i] = FOIntegrate.intOnTriangleRefElement(domainWF.getCompiledLHS()[j][i], 
							new AssembleParam(e, i+1, j+1), params, coords.length, 2);//Laplace Test: 2=80.839 3=80.966, 4=80.967
				}
				b[j] = FOIntegrate.intOnTriangleRefElement(domainWF.getCompiledRHS()[j], 
						new AssembleParam(e, -1, j+1), params, coords.length, 2);
			}
		} else if(domainWF.getFiniteElement().getNumberOfDOFs() == 4) {
			for(int j=0; j<nDOFs; j++) {
				for(int i=0; i<nDOFs; i++) {
					A[j][i] = FOIntegrate.intOnRectangleRefElement(domainWF.getCompiledLHS()[j][i], 
							new AssembleParam(e, i+1, j+1), params, coords.length, 5);
				}
				b[j] = FOIntegrate.intOnRectangleRefElement(domainWF.getCompiledRHS()[j], 
						new AssembleParam(e, -1, j+1), params, coords.length, 5);
			}
		}

		if(null != this.boundaryWF && e.isBorderElement()) {
			ElementList beList = e.getBorderElements();
			for(int n=1;n<=beList.size();n++) {
				//Boundary element
				Element be = beList.at(n);

				//Check node type
				NodeType nodeType = be.getBorderNodeType();
				if(nodeType == NodeType.Neumann || nodeType == NodeType.Robin) {

					double[] beCoords = be.getNodeCoords();
					System.arraycopy(beCoords, 0, beParams, 0, beCoords.length);

					//Update Jacobian on boundary element
					CompiledFunc funcBeJac = this.boundaryWF.getCompiledJac();
					if(null != funcBeJac) 
						funcBeJac.apply(beParams);

					for(int j=0;j<nBeDOFs;j++) {
						for(int i=0;i<nBeDOFs;i++) {
							beA[j][i] = FOIntegrate.intOnLinearRefElement(boundaryWF.getCompiledLHS()[j][i], 
									new AssembleParam(e, i+1, j+1), beParams, beCoords.length, 5);
						}
						beb[j] = FOIntegrate.intOnLinearRefElement(boundaryWF.getCompiledRHS()[j], 
								new AssembleParam(e, -1, j+1), beParams, beCoords.length, 5);
					}
				}
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
		
		FiniteElement fe = this.domainWF.getFiniteElement();
		FiniteElement bfe = fe.getBoundaryFE();
		
		for(Element e : eList) {
			assembleLocal(e);

			for(int j=0;j<nDOFs;j++) {
				int nGlobalRow = fe.getGlobalIndex(mesh, e, j+1);
				for(int i=0;i<nDOFs;i++) {
					int nGlobalCol = fe.getGlobalIndex(mesh, e, i+1);
					stiff.add(nGlobalRow, nGlobalCol, A[j][i]);
				}
				//Local load vector
				load.add(nGlobalRow, b[j]);
			}
			
			if(null != this.boundaryWF && e.isBorderElement()) {
				ElementList beList = e.getBorderElements();
				for(int n=1;n<=beList.size();n++) {
					//Boundary element
					Element be = beList.at(n);

					//Check node type
					NodeType nodeType = be.getBorderNodeType();
					if(nodeType == NodeType.Neumann || nodeType == NodeType.Robin) {
						for(int j=0;j<nBeDOFs;j++) {
							int nGlobalRow = bfe.getGlobalIndex(mesh, be, j+1);
							for(int i=0;i<nBeDOFs;i++) {
								int nGlobalCol = bfe.getGlobalIndex(mesh, be, i+1);
								stiff.add(nGlobalRow, nGlobalCol, beA[j][i]);
							}
							//Local load vector
							load.add(nGlobalRow, beb[j]);
						}
					}
				}
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

	public double[][] getLocalBoundaryStiffMatrix() {
		return beA;
	}
	
	public double[] getLocalBoundaryLoadVector() {
		return beb;
	}

	public Matrix getGlobalStiffMatrix() {
		return gA;
	}

	public Vector getGlobalLoadVector() {
		return gb;
	}
}
