package edu.uta.futureye.test;

import java.util.HashMap;

import edu.uta.futureye.algebra.SolverJBLAS;
import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.application.ModelDOT;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.NodeType;
import edu.uta.futureye.core.intf.Assembler;
import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.basic.FX;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.io.MeshReader;
import edu.uta.futureye.io.MeshWriter;
import edu.uta.futureye.lib.assembler.AssemblerScalar;
import edu.uta.futureye.lib.element.FEBilinearRectangle;
import edu.uta.futureye.lib.weakform.WeakFormLaplace2D;
import edu.uta.futureye.tutorial.Tools;
import edu.uta.futureye.util.container.ElementList;

public class DescreteDerivativeTest {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		MeshReader reader = new MeshReader("rectangle.grd");
		Mesh mesh = reader.read2DMesh();
		mesh.computeNodeBelongsToElements();
		HashMap<NodeType, Function> mapNTF = new HashMap<NodeType, Function>();
		mapNTF.put(NodeType.Dirichlet, null);		

		mesh.markBorderNode(mapNTF);
		//Use element library
		ElementList eList = mesh.getElementList();
		FEBilinearRectangle bilinearRectangle = new FEBilinearRectangle();
		for(int i=1;i<=eList.size();i++)
			bilinearRectangle.assignTo(eList.at(i));
		
		//User defined weak form of PDE (including bounder conditions)
		WeakFormLaplace2D weakForm = new WeakFormLaplace2D();
		
		//-\Delta{u} = f
		//u(x,y)=0, (x,y)\in\partial{\Omega}
		//u=(x^2-9)*(y^2-9)
		//f=-2*(x^2+y^2)+36
		weakForm.setF(FC.c(-2.0).M(
				FX.fx.M(FX.fx).A(FX.fy.M(FX.fy))
			).A(FC.c(36.0))
		);
		
		weakForm.setParam(
				null,
				null,
				new FC(0.05),null //Robin: d*u + k*u_n = q
				); 	
		
		Assembler assembler = new AssemblerScalar(mesh, weakForm);
		assembler.assemble();
		Matrix stiff = assembler.getStiffnessMatrix();
		Vector load = assembler.getLoadVector();
		assembler.imposeDirichletCondition(new FC(0.0));
		
		SolverJBLAS solver = new SolverJBLAS();
		Vector u = solver.solveDGESV(stiff, load);
	    System.out.println("u=");
	    for(int i=1;i<=u.getDim();i++)
	        System.out.println(String.format("%.3f", u.get(i)));
	   
	    MeshWriter writer = new MeshWriter(mesh);
	    writer.writeTechplot("descrete_derivative_test_u.dat", u);
	    
	    Vector ux = Tools.computeDerivative(mesh, u, "x");
	    writer.writeTechplot("descrete_derivative_test_ux.dat", ux);
	    Vector ux2 = edu.uta.futureye.application.Tools.computeDerivative(mesh, u, "x");
	    writer.writeTechplot("descrete_derivative_test_ux2.dat", ux2);
	    
	}

}
