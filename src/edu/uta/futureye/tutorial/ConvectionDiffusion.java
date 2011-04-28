package edu.uta.futureye.tutorial;

import java.util.HashMap;

import edu.uta.futureye.algebra.SolverJBLAS;
import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.NodeType;
import edu.uta.futureye.function.AbstractFunction;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.basic.SpaceVectorFunction;
import edu.uta.futureye.function.basic.Vector2Function;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.function.intf.VectorFunction;
import edu.uta.futureye.io.MeshReader;
import edu.uta.futureye.lib.assembler.AssemblerScalar;
import edu.uta.futureye.lib.element.FELinearTriangle;
import edu.uta.futureye.lib.weakform.WeakFormConvectionDiffusion;
import edu.uta.futureye.util.Constant;
import edu.uta.futureye.util.container.ElementList;

/**
 * Convection-diffusion equation
 *   \frac{\partial{c}}{\partial{t}} = \Nabla{k*\Nabla{c}} - \mathbf{v}\dot\Nabla{c} + f
 * 
 * Time discrete form:
 *   Let:
 *   \frac{\partial{c}}{\partial{t}} = (c_n+1 - c_n)/Dt
 *   We have,
 *   -Dt*\Nabla{k*\Nabla{c_n+1}} + Dt*\mathbf{v}\dot\Nabla{c_n+1} + c_n+1 = Dt*f + c_n
 * 
 * Weak form:
 *   Let c_n+1 := u
 *   Dt*(k*\Nabla{u},\Nabla{w}) + Dt*( (v1*u_x,w)+(v2*u_y,w)+(v3*u_z,w) ) + b*(u,w) = (Dt*f + c_n,w)
 *   
 * where
 *   c=c(x,y,z,t): particles or energy(e.g. salt density, Heat...) are transferred inside 
 *                 a physical system due to two processes: diffusion and convection
 *   k=k(x,y,z): the diffusion coefficient
 *   \mathbf{v}=(v1,v2,v3)':  the convection velocity vector
 *   f=f(x,y,z): the source term
 *   Dt: the time step size
 *   b: b=1 or b=b(x,y,z), if the equation has the term a*c, where a(x,y,z)=b(x,y,z)-1
 * 
 * Boundary condition
 *   c = c0,                  on \Gamma1 (Dirichlet)
 *   d*c + k*c_n = q,         on \Gamma2 (Robin)
 *
 * The following weak form just gives one step computation of c from c_n to c_n+1. 
 *   
 * @author liuyueming
 *
 */
public class ConvectionDiffusion {
	protected String outputFolder = "tutorial\\ConvectionDiffusion";
	protected Mesh mesh = null;
	
	//ConvectionDiffusion weak form
	WeakFormConvectionDiffusion weakForm = new WeakFormConvectionDiffusion();
	
	VectorFunction v = null;
	
	//Time step size
	double Dt;
	
	public void readMesh() {
		//Read a triangle mesh from an input file
		//MeshReader reader = new MeshReader("triangle.grd");
		MeshReader reader = new MeshReader("triangle_refine2.grd");
		mesh = reader.read2DMesh();
		//Geometry relationship
		mesh.computeNodeBelongsToElements();
	}
	
	public void initParam() {
		//Mark border type
		HashMap<NodeType, Function> mapNTF = new HashMap<NodeType, Function>();
		mapNTF.put(NodeType.Dirichlet, new AbstractFunction("x","y"){
			@Override
			public double value(Variable v) {
				//double x = v.get("x");
				double y = v.get("y");
				if(Math.abs(y) < Constant.eps || Math.abs(y-3.0) < Constant.eps)
					return 1;
				else
					return 0;
			}
		});
		mapNTF.put(NodeType.Robin, null);
		mesh.markBorderNode(mapNTF);
		
		//Use element library to assign degree of freedom (DOF) to element
		ElementList eList = mesh.getElementList();
		FELinearTriangle linearTriangle = new FELinearTriangle();
		for(int i=1;i<=eList.size();i++)
			linearTriangle.assignTo(eList.at(i));
	}
	
	public Vector solverOneStep(int step, Function c_n) {
		weakForm.setF(FC.c0);
		weakForm.setParam(FC.c1, FC.c1, c_n, Dt);
		weakForm.setConvectionVelocity(v);
		weakForm.setRobin(FC.c0, FC.c1);
		
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
	    
		Tools.plotVector(mesh, outputFolder, String.format("ConDiff_t%03d.dat",step), u);
	    return u;
		
	}
	
	public void run() {
		readMesh();
		initParam();
		
		//Time step size
		Dt = 0.02;
		
		//the convection velocity vector
		v = new SpaceVectorFunction(2);
		v.set(1, new AbstractFunction("x","y") {
			@Override
			public double value(Variable v) {
				double y = v.get("y");
				return 2*(3.0-Math.abs(y));
			}
		});
		v.set(2, FC.c0);
		
		Tools.plotFunction(mesh, outputFolder, String.format("v1.dat"), v.get(1));
		Tools.plotFunction(mesh, outputFolder, String.format("v2.dat"), v.get(2));
		
		//c0
		Function c_n = new AbstractFunction("x","y") {
			//[-3,3]x[-3,3]			
			@Override
			public double value(Variable v) {
				double x = v.get("x");
				double y = v.get("y");
				if(Math.sqrt((x+2.5)*(x+2.5)+(y-1)*(y-1))<0.3)
					return 40.0;
				if(Math.sqrt((x+2.5)*(x+2.5)+(y+1)*(y+1))<0.3)
					return 40.0;
				else
					return 0.0;
			}
		};
		
		for(int i=1;i<=20;i++) {
			Vector rlt = solverOneStep(i, c_n);
			c_n = new Vector2Function(rlt);
		}
	}
	
	public static void main(String[] args) {
		ConvectionDiffusion cd = new ConvectionDiffusion();
		cd.run();
	}

}
