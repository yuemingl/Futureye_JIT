package edu.uta.futureye.tutorial;

import java.util.HashMap;

import edu.uta.futureye.algebra.SolverJBLAS;
import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.NodeType;
import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.basic.FX;
import edu.uta.futureye.function.basic.Vector2Function;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.io.MeshReader;
import edu.uta.futureye.io.MeshWriter;
import edu.uta.futureye.lib.assembler.AssemblerScalar;
import edu.uta.futureye.lib.element.FELinearTriangle;
import edu.uta.futureye.lib.weakform.WeakFormLaplace2D;
import edu.uta.futureye.util.container.ElementList;

/**
 * d1dt(u) - Laplace(u) = f
 * =>
 *  (u_{n+1}-u_{n})/Dt  - Laplace(u) = f
 * =>
 *  -Dt*Laplace(u) + u_{n+1} = Dt*f + u_{n} 
 *  
 *  u(t=0)=0;
 *  u(x,t)=0, x on border of \Omega
 * 
 * @author liuyueming
 *
 */
public class HeatTransfer_Jiang {
	protected Mesh mesh = null;
	
	//Laplace2D weak form
	WeakFormLaplace2D weakForm = new WeakFormLaplace2D();
	
	//Source term
	Function f = null;
	
	//Time step size
	double Dt;
	
	public void readMesh() {
		//Read a triangle mesh from an input file
		MeshReader reader = new MeshReader("triangle_refine.grd");
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
		HashMap<NodeType, Function> mapNTF = new HashMap<NodeType, Function>();
		mapNTF.put(NodeType.Dirichlet, null);
		mesh.markBorderNode(mapNTF);
		
		//Use element library to assign degree of freedom (DOF) to element
		ElementList eList = mesh.getElementList();
		FELinearTriangle linearTriangle = new FELinearTriangle();
		for(int i=1;i<=eList.size();i++)
			linearTriangle.assignTo(eList.at(i));
		

	}
	
	public Vector solverOneStep(int step, Function u_n) {
		FC FDt = new FC(Dt);
		
		if(step%4==1) {
			f = FC.c(-2.0)
			.M( 
				FX.fx.M(FX.fx).A(FX.fy.M(FX.fy)) )
			.A(
				FC.c(36.0)
			);
//			f = new AbstractFunction("x","y") {
//				@Override
//				public double value(Variable v) {
//					double x = v.get("x");
//					double y = v.get("y");
//					if(Math.abs(x)<Constant.meshEps+0.1 && Math.abs(y)<Constant.meshEps+0.1) {
//						return 36;
//					} else
//						return 0.0;
//				}
//			};
			weakForm.setF(FDt.M(f).A(u_n));
		} else {
			weakForm.setF(u_n);
		}
		
		
		weakForm.setParam(FDt, new FC(1.0), null, null);
		
		//Assemble
		AssemblerScalar assembler = new AssemblerScalar(mesh, weakForm);
		System.out.println("Begin Assemble...");
		assembler.assemble();
		Matrix stiff = assembler.getStiffnessMatrix();
		Vector load = assembler.getLoadVector();
		//Boundary condition
		assembler.imposeDirichletCondition(new FC(0.0));
		System.out.println("Assemble done!");
		
		SolverJBLAS solver = new SolverJBLAS();
		Vector u = solver.solveDGESV(stiff, load);
		System.out.println("u=");
		for(int i=1;i<=u.getDim();i++)
			System.out.println(String.format("%.3f", u.get(i)));	
	    
	    MeshWriter writer = new MeshWriter(mesh);
	    writer.writeTechplot(String.format("tuitoral_HeatTranfer_t%02d.dat",step), u);
	    return u;
		
	}
	
	public void run() {
		readMesh();
		initParam();
		
		//Time step size
		Dt = 0.2;
		
		Function u_n = new FC(0.0);
		for(int i=1;i<=40;i++) {
			Vector rlt = solverOneStep(i, u_n);
			u_n = new Vector2Function(rlt);
		}
	}
	
	public static void main(String[] args) {
		HeatTransfer htf = new HeatTransfer();
		htf.run();
	}
}
