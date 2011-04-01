package edu.uta.futureye.tutorial;

import java.util.HashMap;

import edu.uta.futureye.algebra.Solver;
import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.core.DOF;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.NodeType;
import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.basic.FX;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.io.MeshReader;
import edu.uta.futureye.io.MeshWriter;
import edu.uta.futureye.lib.assembler.AssemblerScalar;
import edu.uta.futureye.lib.shapefun.SFLinearLocal2D;
import edu.uta.futureye.lib.weakform.WeakFormLaplace2D;

/**
 * Problem:
 *   -\Delta{u} = f
 *   u(x,y)=0, (x,y) \in \partial{\Omega}
 * where
 *   \Omega = [-3,3]*[-3,3]
 *   f = -2*(x^2+y^2)+36
 * Solution:
 *   u = (x^2-9)*(y^2-9)
 * 
 * @author Yueming Liu
 */
public class RawLaplace {
	public static void triangle() {
		//Read a triangle mesh from an input file
		MeshReader reader = new MeshReader("triangle.grd");
		Mesh mesh = reader.read2DMesh();
		
		//Geometry relationship
		mesh.computeNodeBelongsToElements();
		
		//Mark border type
		HashMap<NodeType, Function> mapNTF = new HashMap<NodeType, Function>();
		mapNTF.put(NodeType.Dirichlet, null);		
		mesh.markBorderNode(mapNTF);

		//Create 2D linear triangle shape function
		SFLinearLocal2D[] shapeFun = new SFLinearLocal2D[3];
		for(int i=0;i<3;i++)
			shapeFun[i] = new SFLinearLocal2D(i+1);
		
		//Assign degree of freedom(DOF) to element
		for(int i=1;i<=mesh.getElementList().size();i++) {
			Element e = mesh.getElementList().at(i);
			for(int j=1;j<=e.nodes.size();j++) {
				//Create degree of freedom(DOF) object
				DOF dof = new DOF(j,e.nodes.at(j).globalIndex,shapeFun[j-1]);
				e.addNodeDOF(j, dof);
			}
		}
		
		//Laplace2D weak form
		WeakFormLaplace2D weakForm = new WeakFormLaplace2D();
		
		//Right hand side(RHS): f = -2*(x^2+y^2)+36
		Function fx = new FX("x");
		Function fy = new FX("y");
		weakForm.setF(FC.c(-2.0).M( fx.M(fx).A(fy.M(fy)) ).A(FC.c(36.0)));
		
		//Assemble
		AssemblerScalar assembler = new AssemblerScalar(mesh, weakForm);
		System.out.println("Begin Assemble...");
		assembler.assemble();
		Matrix stiff = assembler.getStiffnessMatrix();
		Vector load = assembler.getLoadVector();
		//Boundary condition
		assembler.imposeDirichletCondition(new FC(0.0));
		System.out.println("Assemble done!");
		
		Solver solver = new Solver();
		Vector u = solver.solve(stiff, load);
		System.out.println("u=");
		for(int i=1;i<=u.getDim();i++)
			System.out.println(String.format("%.3f", u.get(i)));	
	    
	    MeshWriter writer = new MeshWriter(mesh);
	    writer.writeTechplot("tuitoral_Laplace.dat", u);
		
	}
	
	public static void main(String[] args) {
		triangle();
	}	
}
