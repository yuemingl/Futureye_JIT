package edu.uta.futureye.tutorial;

import static edu.uta.futureye.function.FMath.C0;
import static edu.uta.futureye.function.FMath.grad;
import static edu.uta.futureye.function.FMath.x;
import static edu.uta.futureye.function.FMath.y;

import java.util.HashMap;

import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.algebra.solver.external.SolverJBLAS;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.NodeType;
import edu.uta.futureye.core.intf.FiniteElement;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.io.MeshReader;
import edu.uta.futureye.io.MeshWriter;
import edu.uta.futureye.lib.assembler.Assembler;
import edu.uta.futureye.lib.element.FELinearTriangle;
import edu.uta.futureye.lib.weakform.WeakForm;
import edu.uta.futureye.util.Utils;

/**
 * Use assembleGlobal() for simple case
 * 
 * <blockquote>
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

public class LaplaceGlobalAssemble {
	public void run() {
		// 1.Read in mesh
		MeshReader reader = new MeshReader("grids/triangle.grd");
		Mesh mesh = reader.read2DMesh();
		// Compute geometry relationship between nodes and elements
		mesh.computeNodeBelongsToElements();

		// 2.Mark boundary type(s)
		HashMap<NodeType, MathFunc> mapNTF = new HashMap<NodeType, MathFunc>();
		mapNTF.put(NodeType.Dirichlet, null); //null => mark all boundary nodes
		mesh.markBorderNode(mapNTF);

		// 3. Weak form definition
		FiniteElement fe = new FELinearTriangle(); // Linear triangular finite element
		final MathFunc f = -2 * (x * x + y * y) + 36; //Right hand side (RHS)
		WeakForm wf = new WeakForm(fe,
				(u,v) -> grad(u, "x", "y").dot(grad(v, "x", "y")), 
				v -> f * v
			);
		wf.compile();

		// 4. Assembly and boundary condition(s)
		Assembler assembler = new Assembler(mesh, wf);
		assembler.assembleGlobal();
		Matrix stiff = assembler.getGlobalStiffMatrix();
		Vector load = assembler.getGlobalLoadVector();
		// Apply zero Dirichlet boundary condition
		Utils.imposeDirichletCondition(stiff, load, fe, mesh, C0);

		// 5. Solve the linear system
		SolverJBLAS solver = new SolverJBLAS();
		Vector u = solver.solveDGESV(stiff, load);
		System.out.println("u=");
		for (int i = 1; i <= u.getDim(); i++)
			System.out.println(String.format("%.3f ", u.get(i)));

		// 6. Output the result to a file with Techplot format
		MeshWriter writer = new MeshWriter(mesh);
		writer.writeTechplot("./tutorial/Laplace2D.dat", u);
	}

	public static void main(String[] args) {
		LaplaceGlobalAssemble ex1 = new LaplaceGlobalAssemble();
		ex1.run();
	}
}
