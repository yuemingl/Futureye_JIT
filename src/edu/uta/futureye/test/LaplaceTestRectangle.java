package edu.uta.futureye.test;

import static edu.uta.futureye.function.FMath.C0;
import static edu.uta.futureye.function.FMath.grad;
import static edu.uta.futureye.function.FMath.x;
import static edu.uta.futureye.function.FMath.y;

import java.util.HashMap;

import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.algebra.solver.external.SolverJBLAS;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.NodeType;
import edu.uta.futureye.core.intf.FiniteElement;
import edu.uta.futureye.function.MultiVarFunc;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.io.MeshReader;
import edu.uta.futureye.io.MeshWriter;
import edu.uta.futureye.lib.assembler.Assembler;
import edu.uta.futureye.lib.assembler.BasicAssembler;
import edu.uta.futureye.lib.assembler.DomainBoundaryAssemblerRaw;
import edu.uta.futureye.lib.element.FEBilinearRectangle;
import edu.uta.futureye.lib.weakform.WeakForm;
import edu.uta.futureye.util.Utils;

/**
 * <blockquote><pre>
 * Solve
 *   -k*\Delta{u} = f  in \Omega
 *   u(x,y) = 0,       on boundary x=3.0 of \Omega
 *   k*u_n + u = 0.01,   on other boundary of \Omega
 * where
 *   \Omega = [-3,3]*[-3,3]
 *   k = 2
 *   f = -4*(x^2+y^2)+72
 *   u_n = \frac{\partial{u}}{\partial{n}}
 *   n: outer unit normal vector of \Omega
 * </blockquote></pre>
 */

public class LaplaceTestRectangle {
	public void run() {
		// 1.Read mesh
		MeshReader reader = new MeshReader("grids/rectangle.grd");
		Mesh mesh = reader.read2DMesh();
		// Compute geometry relationship between nodes and elements
		mesh.computeNodeBelongsToElements();

		// 2.Mark border types
		HashMap<NodeType, MathFunc> mapNTF = new HashMap<NodeType, MathFunc>();
		//Robin type on boundary x=3.0 of \Omega
		mapNTF.put(NodeType.Robin, new MultiVarFunc("Robin", "x","y"){
			@Override
			public double apply(double... args) {
				if(Math.abs(3.0-args[this.argIdx[0]]) < 0.01)
					return 1.0; //this is Robin condition
				else
					return -1.0;
			}
		});
		//Dirichlet type on other boundary of \Omega
		mapNTF.put(NodeType.Dirichlet, null);
		mesh.markBorderNode(mapNTF);

		// 3.Use finite element library to assign degrees of
		// freedom (DOF) to element
		FEBilinearRectangle fe = new FEBilinearRectangle();
//		for(Element e : mesh.getElementList())
//			fet.assignTo(e);

		//4. Weak forms
		//Right hand side(RHS):
		final MathFunc f = - 4*(x*x + y*y) + 72;
		//Weak form in the domain
		WeakForm wf = new WeakForm(fe, 
				(u,v) -> 2*grad(u, "x", "y").dot(grad(v, "x", "y")), 
				v -> f * v
			);
		wf.compile();
		//Weak form on the boundary (robin condition)
		WeakForm bwf = new WeakForm(
				fe.getBoundaryFE(), //boundary finite element
				(u,v) -> 2*u*v, 
				v -> 0.01 * v
			);
		bwf.compile();

		// 5.Assembly process
		//Assembler assembler = new Assembler(wf, bwf);
		DomainBoundaryAssemblerRaw assembler = new DomainBoundaryAssemblerRaw(wf, bwf);
		assembler.assembleGlobal(mesh);
		Matrix stiff = assembler.getGlobalStiffMatrix();
		Vector load = assembler.getGlobalLoadVector();
		// Boundary condition
		Utils.imposeDirichletCondition(stiff, load, fe, mesh, C0);

		// 6.Solve linear system
		SolverJBLAS solver = new SolverJBLAS();
		Vector u = solver.solveDGESV(stiff, load);
		System.out.println("u=");
		for (int i = 1; i <= u.getDim(); i++)
			System.out.println(String.format("%.3f ", u.get(i)));

		// 7.Output results to an Techplot format file
		MeshWriter writer = new MeshWriter(mesh);
		writer.writeTechplot("Laplace2DRectangle.dat", u);
	}

	
	/**
	 * Use two basic assemblers
	 */
	public void run2() {
		// 1.Read mesh
		MeshReader reader = new MeshReader("grids/rectangle.grd");
		Mesh mesh = reader.read2DMesh();
		// Compute geometry relationship between nodes and elements
		mesh.computeNodeBelongsToElements();

		// 2.Mark border types
		HashMap<NodeType, MathFunc> mapNTF = new HashMap<NodeType, MathFunc>();
		//Robin type on boundary x=3.0 of \Omega
		mapNTF.put(NodeType.Robin, new MultiVarFunc("Robin", "x","y"){
			@Override
			public double apply(double... args) {
				if(Math.abs(3.0-args[this.argIdx[0]]) < 0.01)
					return 1.0; //this is Robin condition
				else
					return -1.0;
			}
		});
		//Dirichlet type on other boundary of \Omega
		mapNTF.put(NodeType.Dirichlet, null);
		mesh.markBorderNode(mapNTF);

		// 3.Use finite element library to assign degrees of
		// freedom (DOF) to element
		FEBilinearRectangle fe = new FEBilinearRectangle();
//		for(Element e : mesh.getElementList())
//			fet.assignTo(e);

		//4. Weak forms
		//Right hand side(RHS):
		final MathFunc f = - 4*(x*x + y*y) + 72;
		//Weak form in the domain
		WeakForm wf = new WeakForm(fe, 
				(u,v) -> 2*grad(u, "x", "y").dot(grad(v, "x", "y")), 
				v -> f * v
			);
		wf.compile();
		//Weak form on the boundary (robin condition)
		FiniteElement beFE = fe.getBoundaryFE();
		WeakForm bwf = new WeakForm(beFE, //boundary finite element
				(u,v) -> 2*u*v, 
				v -> 0.01 * v
			);
		bwf.compile();

		// 5.Assembly process
		BasicAssembler domainAssembler = new BasicAssembler(mesh, wf);
		// Assemble on the elements
		domainAssembler.assembleGlobal();
		Matrix stiff = domainAssembler.getGlobalStiffMatrix();
		Vector load = domainAssembler.getGlobalLoadVector();
		// Assemble on the boundary elements
		BasicAssembler boundaryAssembler = new BasicAssembler(mesh, bwf);
		int beNumDOFs = fe.getBoundaryFE().getNumberOfDOFs();
		for (Element e : mesh.getElementList()) {
			// Use BasicAssembler to assemble boundary elements
			for(Element be : e.getBorderElements()) {
				// Check node type
				NodeType nodeType = be.getBorderNodeType();
				if(nodeType == NodeType.Neumann || nodeType == NodeType.Robin) {
					// Assemble locally on boundary element
					boundaryAssembler.assembleLocal(be);
					double[][] beA = boundaryAssembler.getLocalStiffMatrix();
					double[] beb = boundaryAssembler.getLocalLoadVector();
					
					// Get all local DOF of the boundary element for local-global indexing
					for (int j = 0; j < beNumDOFs; j++) {
						int nGlobalRow = beFE.getGlobalIndex(mesh, be, j+1);
						for (int i = 0; i < beNumDOFs; i++) {
							int nGlobalCol = beFE.getGlobalIndex(mesh, be, i+1);
							stiff.add(nGlobalRow, nGlobalCol, beA[j][i]);
						}
						load.add(nGlobalRow, beb[j]);
					}
				}
			}
		}
		// Boundary condition
		Utils.imposeDirichletCondition(stiff, load, fe, mesh, C0);

		// 6.Solve linear system
		SolverJBLAS solver = new SolverJBLAS();
		Vector u = solver.solveDGESV(stiff, load);
		System.out.println("u=");
		for (int i = 1; i <= u.getDim(); i++)
			System.out.println(String.format("%.3f ", u.get(i)));

		// 7.Output results to an Techplot format file
		MeshWriter writer = new MeshWriter(mesh);
		writer.writeTechplot("Laplace2DRectangle.dat", u);
	}
	
	
	/**
	 * Use two chained assemblers
	 */
	public void run3() {
		// 1.Read mesh
		MeshReader reader = new MeshReader("grids/rectangle.grd");
		Mesh mesh = reader.read2DMesh();
		// Compute geometry relationship between nodes and elements
		mesh.computeNodeBelongsToElements();

		// 2.Mark border types
		HashMap<NodeType, MathFunc> mapNTF = new HashMap<NodeType, MathFunc>();
		//Robin type on boundary x=3.0 of \Omega
		mapNTF.put(NodeType.Robin, new MultiVarFunc("Robin", "x","y"){
			@Override
			public double apply(double... args) {
				if(Math.abs(3.0-args[this.argIdx[0]]) < 0.01)
					return 1.0; //this is Robin condition
				else
					return -1.0;
			}
		});
		//Dirichlet type on other boundary of \Omega
		mapNTF.put(NodeType.Dirichlet, null);
		mesh.markBorderNode(mapNTF);

		// 3.Use finite element library to assign degrees of
		// freedom (DOF) to element
		FEBilinearRectangle fe = new FEBilinearRectangle();
//		for(Element e : mesh.getElementList())
//			fet.assignTo(e);

		//4. Weak forms
		//Right hand side(RHS):
		final MathFunc f = - 4*(x*x + y*y) + 72;
		//Weak form in the domain
		WeakForm wf = new WeakForm(fe, 
				(u,v) -> 2*grad(u, "x", "y").dot(grad(v, "x", "y")), 
				v -> f * v
			);
		wf.compile();
		//Weak form on the boundary (robin condition)
		WeakForm bwf = new WeakForm(
				fe.getBoundaryFE(), //boundary finite element
				(u,v) -> 2*u*v, 
				v -> 0.01 * v
			);
		bwf.compile();

		// 5.Assembly process
		Assembler domainAssembler = new Assembler(mesh, wf);
		domainAssembler.assembleGlobal();
		Assembler boundaryAssembler = new Assembler(domainAssembler, bwf);
		for (Element e : mesh.getElementList()) {
			// Use BasicAssembler to assemble boundary elements
			for(Element be : e.getBorderElements()) {
				// Check node type
				NodeType nodeType = be.getBorderNodeType();
				if(nodeType == NodeType.Neumann || nodeType == NodeType.Robin) {
					// Assemble boundary element into global stiff matrix and load vector
					boundaryAssembler.assembleGlobal(be);
				}
			}
		}
		Matrix stiff = boundaryAssembler.getGlobalStiffMatrix();
		Vector load = boundaryAssembler.getGlobalLoadVector();
		// Boundary condition
		Utils.imposeDirichletCondition(stiff, load, fe, mesh, C0);

		// 6.Solve linear system
		SolverJBLAS solver = new SolverJBLAS();
		Vector u = solver.solveDGESV(stiff, load);
		System.out.println("u=");
		for (int i = 1; i <= u.getDim(); i++)
			System.out.println(String.format("%.3f ", u.get(i)));

		// 7.Output results to an Techplot format file
		MeshWriter writer = new MeshWriter(mesh);
		writer.writeTechplot("Laplace2DRectangle.dat", u);
	}
	public static void main(String[] args) {
		LaplaceTestRectangle ex1 = new LaplaceTestRectangle();
		ex1.run();  //23.518
		ex1.run2(); //23.518
		ex1.run3(); //23.518
	}
}
