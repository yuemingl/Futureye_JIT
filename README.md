FuturEye_JIT is a Finite Element Method (FEM) library for solving Partial Defferential Equation (PDE) based forward/inverse problems. Java is widely used in the industry level software and systems. The JVM's HotSpot Just-In-Time (JIT) compiler allows the speed of the Java program to approach that of a native application. FuturEye_JIT makes use of the JIT feature by generating small sized Java bytecode functions at runtime for the mathematical expressions in the definitoin of a problem. The resulting code is as efficient as hand written FORTRAN FEM code.

FuturEye_JIT provides a mathematically appealing way for building functions and weak forms of a problem by using the Java operator overloading technique (https://github.com/amelentev/java-oo). Thus a concise, natrual and human friendly application programming interface (API) is provided. 

The basic components in FEM are abstracted out, such as node, element, mesh, degree of freedom and shape function. The data structure and operation of these classes are encapsulated properly by OOP. The things that make FuturEye_JIT unique among the existing object-oriented FEM software or libraries are the function classes. The behavior of the function classes in Futureye_JIT is very similar to that in the human readable mathematical context. For example algebra of functions, function derivatives and composition of functions. Especially in FEM theory, shape functions, Jacobin of coordinate transform and numerical integration are all based on the function classes. This feature leads to a close integration between the mathematical theory and it's computer implementation. FuturEye_JIT is designed to solve 1D, 2D and 3D PED problems with scalar or vector valued unknowns. It is motivated by solving PDE based inverse problem in the application of optical tomography where the word 'FuturEye' came from.

FuturEye_JIT is suitable for various purposes:

*Teaching: The feature of close relation to the mathematical theory of FEM helps the students to understand basic FEM concepts, e.g. shape functions, Jacobian and the assembly process.

*Research: FuturEye_JIT helps researchers quickly develop and test their models, experiment with data and algorithms. e.g. new equations, new finite elements and new solution methods without concerning too much about programming basic components in FEM.

*Engineering: Futureye_JIT is designed to be robust and efficient. Industry level of applications can be easily built with it. 

### Laplace Example ###
<img src='https://lh6.googleusercontent.com/_Cil2MFH7iLM/TN19jeWDEdI/AAAAAAAAABg/WI64bT_jUAY/s800/FutureEyeFirstTest2.png.jpg' />

### Solution ###
Mesh, solution contour and 3D view:
<img src='https://lh5.googleusercontent.com/_Cil2MFH7iLM/TN19jH3fdUI/AAAAAAAAABc/bjKllifWW0g/s288/FutureEyeFirstTest.png.jpg' /><img src='https://lh3.googleusercontent.com/_Cil2MFH7iLM/TN19j0Dy4pI/AAAAAAAAABk/OTdlyX_Paio/s288/FutureEyeFirstTest3D.png.jpg' />

### Code ###

```java
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
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.io.MeshReader;
import edu.uta.futureye.io.MeshWriter;
import edu.uta.futureye.lib.assembler.Assembler;
import edu.uta.futureye.lib.element.FELinearTriangle;
import edu.uta.futureye.lib.weakform.WeakForm;
import edu.uta.futureye.util.Utils;

public class LaplaceGlobalAssemble {
	public void run() {
		// 1.Read in mesh
		MeshReader reader = new MeshReader("grids/triangle.grd");
		Mesh mesh = reader.read2DMesh();
		// Compute geometry relationship between nodes and elements
		mesh.computeNodeBelongsToElements();

		// 2.Mark boundary type(s)
		HashMap<NodeType, MathFunc> mapNTF = new HashMap<NodeType, MathFunc>();
		mapNTF.put(NodeType.Dirichlet, null);
		mesh.markBorderNode(mapNTF);

		// 3. Weak form definition
		final MathFunc f = -2 * (x * x + y * y) + 36; //Right hand side (RHS)
		WeakForm wf = new WeakForm(
				new FELinearTriangle(), // Linear triangular finite element
				(u,v) -> grad(u, "x", "y").dot(grad(v, "x", "y")), 
				v -> f * v
			);
		wf.compile();

		// 4. Assembly and boundary condition(s)
		Assembler assembler = new Assembler(wf);
		assembler.assembleGlobal(mesh);
		Matrix stiff = assembler.getGlobalStiffMatrix();
		Vector load = assembler.getGlobalLoadVector();
		// Apply zero Dirichlet boundary condition
		Utils.imposeDirichletCondition(stiff, load, mesh, C0);

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
```
