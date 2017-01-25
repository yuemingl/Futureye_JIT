package edu.uta.futureye.lib.assembler;

import edu.uta.futureye.algebra.SparseMatrixRowMajor;
import edu.uta.futureye.algebra.SparseVectorHashMap;
import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.core.DOF;
import edu.uta.futureye.core.DOFOrder;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.NodeType;
import edu.uta.futureye.lib.weakform.WeakForm;
import edu.uta.futureye.util.container.DOFList;
import edu.uta.futureye.util.container.ElementList;

public class Assembler {
	BasicAssembler domainAss;
	BasicAssembler boundaryAss;
	
	Matrix gA;
	Vector gb;
	
	/**
	 * 
	 * @param domainWeakForm
	 * @param boundaryWeakForm
	 */
	public Assembler(WeakForm domainWeakForm, WeakForm boundaryWeakForm) {
		this.domainAss = new BasicAssembler(domainWeakForm);
		this.boundaryAss = new BasicAssembler(boundaryWeakForm);
		
	}
	
	/**
	 * Assemble local stiff and load on a give element
	 * @param e
	 */
	public void assembleLocal(Element e) {
		// Assemble on domain element
		domainAss.assembleLocal(e);
		
		// Assemble on boundary element
		for(Element be : e.getBorderElements()) {
			//Check node type
			NodeType nodeType = be.getBorderNodeType();
			if(nodeType == NodeType.Neumann || nodeType == NodeType.Robin) {
				//Associate boundary FiniteElement object to the boundary element
				this.boundaryAss.weakForm.getFiniteElement().assignTo(be);
					this.boundaryAss.assembleLocal(be);
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
			
			// Assemble locally
			assembleLocal(e);
			
			// Get local-global indexing
			DOFList DOFs = e.getAllDOFList(DOFOrder.NEFV);
			for(int j=0;j<DOFs.size();j++) {
				DOF dofJ = DOFs.at(j+1);
				for(int i=0;i<DOFs.size();i++) {
					DOF dofI = DOFs.at(i+1);
					stiff.add(dofJ.getGlobalIndex(), dofI.getGlobalIndex(), this.domainAss.A[j][i]);
				}
				load.add(dofJ.getGlobalIndex(), this.domainAss.b[j]);
			}
			
			// Use BasicAssembler to assemble boundary elements
			for(Element be : e.getBorderElements()) {
				//Check node type
				NodeType nodeType = be.getBorderNodeType();
				if(nodeType == NodeType.Neumann || nodeType == NodeType.Robin) {
					// Get local-global indexing
					DOFList beDOFs = be.getAllDOFList(DOFOrder.NEFV);
					for(int j=0;j<beDOFs.size();j++) {
						DOF beDOFJ = beDOFs.at(j+1);
						for(int i=0;i<beDOFs.size();i++) {
							DOF beDOFI = beDOFs.at(i+1);
							stiff.add(beDOFJ.getGlobalIndex(), beDOFI.getGlobalIndex(), this.boundaryAss.A[j][i]);
						}
						load.add(beDOFJ.getGlobalIndex(), this.boundaryAss.b[j]);
					}
				}
			}
		}
		//update gA and gb
		this.gA = stiff;
		this.gb = load;
	}
	
	public double[][] getLocalStiffMatrix() {
		return this.domainAss.getLocalStiffMatrix();
	}
	
	public double[] getLocalLoadVector() {
		return this.domainAss.getLocalLoadVector();
	}

	public double[][] getLocalBoundaryStiffMatrix() {
		return this.boundaryAss.getLocalStiffMatrix();
	}
	
	public double[] getLocalBoundaryLoadVector() {
		return this.boundaryAss.getLocalLoadVector();
	}

	public Matrix getGlobalStiffMatrix() {
		return gA;
	}

	public Vector getGlobalLoadVector() {
		return gb;
	}
}
