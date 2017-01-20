package edu.uta.futureye.test;

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
import edu.uta.futureye.core.intf.LHSExpr;
import edu.uta.futureye.core.intf.RHSExpr;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.io.MeshReader;
import edu.uta.futureye.io.MeshWriter;
import edu.uta.futureye.lib.element.AssemblerJIT;
import edu.uta.futureye.lib.element.FELinearTriangle;
import edu.uta.futureye.lib.element.FELinearTriangleJIT;
import edu.uta.futureye.lib.element.WeakFormJIT;
import edu.uta.futureye.util.MeshGenerator;
import edu.uta.futureye.util.Utils;
import edu.uta.futureye.util.container.ElementList;


/**
 * Use assembleGlobal() for simple case
 * 
 * <blockquote><pre>
 * Problem:
 *   -\Delta{u} = f
 *   u(x,y)=0, (x,y) \in \partial{\Omega}
 * where
 *   \Omega = [-3,3]*[-3,3]
 *   f = 2 * pi * pi * sin ( pi * x ) * sin ( pi * y )
 * Solution:
 *   u = (x^2-9)*(y^2-9)
 * </blockquote></pre>
 * 
 * @author liuyueming
 */

public class LaplaceTestJITForPaperRefactoredAssembeGlobal {
	public Mesh mesh; //mesh object
	public Vector u; //solution vector
	
	public void run(int nNodes) {
		int n = 51;
		boolean solveSystem = false;

		// 1.Generate mesh
		Mesh mesh = null;
		if (solveSystem) {
			MeshReader reader = new MeshReader("triangle.grd");
			mesh = reader.read2DMesh();
		} else {
			mesh = MeshGenerator.rectangle(-3, 3, -3, 3, n, n);
		}

		// Compute geometry relationship between nodes and elements
		mesh.computeNodeBelongsToElements();

		// 2.Mark border types
		HashMap<NodeType, MathFunc> mapNTF = new HashMap<NodeType, MathFunc>();
		mapNTF.put(NodeType.Dirichlet, null);
		mesh.markBorderNode(mapNTF);

		// 3.Use element library to assign degrees of
		// freedom (DOF) to element
		ElementList eList = mesh.getElementList();
		FELinearTriangle feLT = new FELinearTriangle();
		for (int i = 1; i <= eList.size(); i++)
			feLT.assignTo(eList.at(i));

		// Finite element
		FELinearTriangleJIT fet = new FELinearTriangleJIT();

		// Right hand side(RHS):
		// MathFunc f = 2.0 * PI * PI * sin ( PI * x ) * sin ( PI * y );
		final MathFunc f = -2 * (x * x + y * y) + 36;

		WeakFormJIT wf = new WeakFormJIT(fet, new LHSExpr() {
			public MathFunc apply(MathFunc u, MathFunc v) {
				return grad(u, "x", "y").dot(grad(v, "x", "y"));
			}
		}, new RHSExpr() {
			public MathFunc apply(MathFunc v) {
				return f * v;
			}
		});

		long startCompile = System.currentTimeMillis();
		wf.compile();
		System.out.println("Compile time: "
				+ (System.currentTimeMillis() - startCompile));

		// 5.Assembly process
		AssemblerJIT assembler = new AssemblerJIT(wf);
		Matrix stiff = null;
		Vector load = null;

		long start = System.currentTimeMillis();
		int NN = nNodes / ((n - 1) * (n - 1)); // 10000*512/eList.size();
		if (solveSystem)
			NN = 1;
		for (int ii = 0; ii < NN; ii++) {
			assembler.assembleGlobal(mesh);
			stiff = assembler.getGlobalStiffMatrix();
			load = assembler.getGlobalLoadVector();
		}
		System.out.println("Nodes=" + nNodes + ", Aassembly time: "
				+ (System.currentTimeMillis() - start) + "ms");

		if (solveSystem) {
			// Boundary condition
			Utils.imposeDirichletCondition(stiff, load, mesh, C0);

			// 6.Solve linear system
			SolverJBLAS solver = new SolverJBLAS();
			Vector u = solver.solveDGESV(stiff, load);
			System.out.println("u=");
			for (int i = 1; i <= u.getDim(); i++)
				System.out.println(String.format("%.3f ", u.get(i)));

			// 7.Output results to an Techplot format file
			MeshWriter writer = new MeshWriter(mesh);
			writer.writeTechplot("./tutorial/Laplace2D.dat", u);

			this.mesh = mesh;
			this.u = u;
		}
	}
	
	public static void main(String[] args) {
		LaplaceTestJITForPaperRefactoredAssembeGlobal ex1 = new LaplaceTestJITForPaperRefactoredAssembeGlobal();
		ex1.run(10000);
		ex1.run(100000);
		ex1.run(10000 * 512);
		ex1.run(1000000);
		ex1.run(10000000);
		/**
		 * 
compileWithASM, assemble globally
Nodes=10000, Aassembly time: 152ms
Nodes=100000, Aassembly time: 481ms
Nodes=5120000, Aassembly time: 23915ms
Nodes=1000000, Aassembly time: 4690ms
Nodes=10000000, Aassembly time: 46628ms

		 */
	}
}
