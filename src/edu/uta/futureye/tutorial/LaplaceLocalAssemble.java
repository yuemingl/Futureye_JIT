package edu.uta.futureye.tutorial;

import static edu.uta.futureye.function.FMath.C0;
import static edu.uta.futureye.function.FMath.grad;
import static edu.uta.futureye.function.FMath.x;
import static edu.uta.futureye.function.FMath.y;

import java.util.HashMap;

import edu.uta.futureye.algebra.SparseMatrixRowMajor;
import edu.uta.futureye.algebra.SparseVectorHashMap;
import edu.uta.futureye.algebra.intf.SparseMatrix;
import edu.uta.futureye.algebra.intf.SparseVector;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.algebra.solver.external.SolverJBLAS;
import edu.uta.futureye.core.DOF;
import edu.uta.futureye.core.DOFOrder;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.NodeType;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.io.MeshReader;
import edu.uta.futureye.io.MeshWriter;
import edu.uta.futureye.lib.assembler.Assembler;
import edu.uta.futureye.lib.element.FELinearTriangle;
import edu.uta.futureye.lib.weakform.WeakForm;
import edu.uta.futureye.util.Utils;
import edu.uta.futureye.util.container.DOFList;

/**
 * Use assembleLocal() to get local stiff matrix and local load vector in an
 * element This gives user the ability to assemble their own global stiff matrix
 * and global load vector <blockquote>
 * 
 * <pre>
 * Problem:
 *   -\Delta{u} = f
 *   u(x,y)=0, (x,y) \in \partial{\Omega}
 * where
 *   \Omega = [-3,3]*[-3,3]
 *   f = -2*(x^2+y^2)+36
 * Solution:
 *   u = (x^2-9)*(y^2-9)
 * </blockquote>
 * </pre>
 * 
 * @author liuyueming
 */

public class LaplaceLocalAssemble {
	public void run() {
		// 1. Read in mesh
		MeshReader reader = new MeshReader("grids/triangle.grd");
		Mesh mesh = reader.read2DMesh();
		// Compute geometry relationship between nodes and elements
		mesh.computeNodeBelongsToElements();

		// 2. Mark boundary types
		HashMap<NodeType, MathFunc> mapNTF = new HashMap<NodeType, MathFunc>();
		mapNTF.put(NodeType.Dirichlet, null); //null => mark all boundary nodes
		mesh.markBorderNode(mapNTF);

		// 3. Weak form definition
		FELinearTriangle fe = new FELinearTriangle(); //Linear triangular finite element
		final MathFunc f = -2 * (x * x + y * y) + 36; //Right hand side (RHS)
		WeakForm wf = new WeakForm(fe, 
				(u,v) -> grad(u, "x", "y").dot(grad(v, "x", "y")), 
				v -> f * v
			);
		wf.compile();

		// 5. Assembly and boundary condition(s)
		Assembler assembler = new Assembler(mesh, wf);
		int dim = mesh.getNodeList().size();
		SparseMatrix stiff = new SparseMatrixRowMajor(dim, dim);
		SparseVector load = new SparseVectorHashMap(dim);
		int nDOFs = fe.getNumberOfDOFs();
		for (Element e : mesh.getElementList()) {
			assembler.assembleLocal(e);
			double[][] A = assembler.getLocalStiffMatrix();
			double[] b = assembler.getLocalLoadVector();
			for (int j = 0; j < nDOFs; j++) {
				int nGlobalRow = fe.getGlobalIndex(mesh, e, j+1);
				for (int i = 0; i < nDOFs; i++) {
					int nGlobalCol = fe.getGlobalIndex(mesh, e, i+1);
					stiff.add(nGlobalRow, nGlobalCol, A[j][i]);
				}
				// Local load vector
				load.add(nGlobalRow, b[j]);
			}
		}
		// Apply zero Dirichlet boundary condition
		Utils.imposeDirichletCondition(stiff, load, fe, mesh, C0);

		// 6. Solve linear system
		SolverJBLAS solver = new SolverJBLAS();
		Vector u = solver.solveDGESV(stiff, load);
		System.out.println("u=");
		for (int i = 1; i <= u.getDim(); i++)
			System.out.println(String.format("%.3f ", u.get(i)));

		// 7. Output the result to a file with Techplot format
		MeshWriter writer = new MeshWriter(mesh);
		writer.writeTechplot("./tutorial/Laplace2D.dat", u);
	}

	public static void main(String[] args) {
		LaplaceLocalAssemble ex1 = new LaplaceLocalAssemble();
		ex1.run();
	}
}
