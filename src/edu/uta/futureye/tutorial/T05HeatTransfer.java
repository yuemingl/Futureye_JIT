package edu.uta.futureye.tutorial;

import java.io.File;
import java.util.HashMap;

import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.algebra.solver.external.SolverJBLAS;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.NodeType;
import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.basic.FX;
import edu.uta.futureye.function.basic.Vector2Function;
import edu.uta.futureye.function.intf.MathFun;
import edu.uta.futureye.io.MeshReader;
import edu.uta.futureye.io.MeshWriter;
import edu.uta.futureye.lib.assembler.AssemblerScalar;
import edu.uta.futureye.lib.element.FELinearTriangle;
import edu.uta.futureye.lib.weakform.WeakFormLaplace2D;
import edu.uta.futureye.util.container.ElementList;

/**
 * <blockquote><pre>
 * Heat transfer problem:
 * d(u)/dt - Laplace(u) = f
 * =>
 *  (u_{n+1}-u_{n})/Dt  - Laplace(u) = f
 * =>
 *  -Dt*Laplace(u) + u_{n+1} = Dt*f + u_{n} 
 *  
 *  u(t=0)=0;
 *  u(x,t)=0, x on border of \Omega
 * <blockquote><pre>
 * 
 * @author liuyueming
 */
public class T05HeatTransfer {
	String outputFolder = "./tutorial/HeatTranfer";
	Mesh mesh = null;
	
	//Laplace2D weak form
	WeakFormLaplace2D weakForm = new WeakFormLaplace2D();
	
	//Source term
	MathFun f = null;
	
	//Time step size
	double Dt;
	
	public void readMesh() {
		//Read a triangle mesh from an input file
		MeshReader reader = new MeshReader("triangle.grd");
		mesh = reader.read2DMesh();
		//Geometry relationship
		mesh.computeNodeBelongsToElements();
	}
	
	public void initParam() {
		//Right hand side(RHS): f = -2*(x^2+y^2)+36
		f = FC.c(-2.0)
			.M( 
				FX.fx.M(FX.fx).A(FX.fy.M(FX.fy)) )
			.A(
				FC.c(36.0)
			);
				
		//Mark border type
		HashMap<NodeType, MathFun> mapNTF = new HashMap<NodeType, MathFun>();
		mapNTF.put(NodeType.Dirichlet, null);
		mesh.markBorderNode(mapNTF);
		
		//Use element library to assign degree of freedom (DOF) to element
		ElementList eList = mesh.getElementList();
		FELinearTriangle linearTriangle = new FELinearTriangle();
		for(int i=1;i<=eList.size();i++)
			linearTriangle.assignTo(eList.at(i));
		
	    File file = new File(outputFolder);
		if(!file.exists()) {
			file.mkdirs();
		}
	}
	
	public Vector solverOneStep(int step, MathFun u_n) {
		FC FDt = new FC(Dt);
		
		weakForm.setF(FDt.M(f).A(u_n));
		weakForm.setParam(FDt, new FC(1.0), null, null);
		
		//Assemble
		AssemblerScalar assembler = new AssemblerScalar(mesh, weakForm);
		assembler.assemble();
		Matrix stiff = assembler.getStiffnessMatrix();
		Vector load = assembler.getLoadVector();
		//Boundary condition
		assembler.imposeDirichletCondition(new FC(0.0));
		
		SolverJBLAS solver = new SolverJBLAS();
		Vector u = solver.solveDGESV(stiff, load);
		System.out.println("u=");
		for(int i=1;i<=u.getDim();i++)
			System.out.println(String.format("%.3f", u.get(i)));	
	    
	    MeshWriter writer = new MeshWriter(mesh);
	    writer.writeTechplot(String.format("%s/u_t%02d.dat",outputFolder, step), u);
	    return u;
		
	}
	
	public void run() {
		readMesh();
		initParam();
		
		//Time step size
		Dt = 0.2;
		
		MathFun u_n = new FC(0.0);
		for(int i=1;i<=25;i++) {
			Vector rlt = solverOneStep(i, u_n);
			u_n = new Vector2Function(rlt);
		}
	}
	
	public static void main(String[] args) {
		T05HeatTransfer htf = new T05HeatTransfer();
		htf.run();
	}
}
