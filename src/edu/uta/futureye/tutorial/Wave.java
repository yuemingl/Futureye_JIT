package edu.uta.futureye.tutorial;

import java.io.File;
import java.util.HashMap;

import edu.uta.futureye.algebra.CompressedRowMatrix;
import edu.uta.futureye.algebra.FullVector;
import edu.uta.futureye.algebra.Solver;
import edu.uta.futureye.algebra.SparseMatrix;
import edu.uta.futureye.algebra.SparseVector;
import edu.uta.futureye.algebra.intf.AlgebraMatrix;
import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.core.NodeType;
import edu.uta.futureye.function.AbstractFunction;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.basic.Vector2Function;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.io.MeshReader;
import edu.uta.futureye.io.MeshWriter;
import edu.uta.futureye.lib.assembler.AssemblerScalar;
import edu.uta.futureye.lib.element.FEBilinearRectangle;
import edu.uta.futureye.lib.element.FELinearTriangle;
import edu.uta.futureye.lib.weakform.WeakFormLaplace2D;
import edu.uta.futureye.util.list.ElementList;
import edu.uta.futureye.util.list.NodeList;

/**
 * d2dt(u) - c2*Laplace(u) = 0
 * => 
 * (u_{n+2}-2*u_{n+1}+u_{n})/Dt^2  - c2*Laplace(u) = 0
 * =>
 * -Dt^2*c2*Laplace(u_{n+2}) + u_{n+2} = 2*u_{n+1} - u_{n} 
 *  
 *  u(t=0) = u0(x);
 *  u + u_n(x,t) = 0, x on border of \Omega
 * 
 * @author liuyueming
 *
 */
public class Wave {
	protected static String outputFolder = "tutorial\\wave";
	protected Mesh mesh = null;
	
	//Laplace2D weak form
	WeakFormLaplace2D weakForm = new WeakFormLaplace2D();
	
	//Time step size
	double Dt;
	
	//Speed^2
	double c2;

	//2*u_{n+1} - u_{n}
	Function u_n = null;
	Function u_n1 = null;
	
	//u(t=0)=u0(x)
	Function u0 = null;
	
	//Function f = null;
	
	public void readMesh() {
		//Read a triangle mesh from an input file
		//[-3,3]x[-3,3]
		//MeshReader reader = new MeshReader("triangle_refine80x80.grd");
		
		//[-10,10]x[-100,10]
		//MeshReader reader = new MeshReader("rectangle_big.grd");
		//MeshReader reader = new MeshReader("rectangle_big2.grd");
		//MeshReader reader = new MeshReader("triangle_big.grd");
		//MeshReader reader = new MeshReader("triangle_big2.grd");
		
		//凸型区域[-40,40]x[-30,25]	
		//MeshReader reader = new MeshReader("two_rectangle.grd");
		MeshReader reader = new MeshReader("three_rectangle.grd");
		
		//[-100,100]x[-100,100]
		//MeshReader reader = new MeshReader("triangle_bigbig.grd");
		//MeshReader reader = new MeshReader("triangle_bigbig80x80.grd");
		mesh = reader.read2DMesh();
		//Geometry relationship
		mesh.computeNodesBelongToElement();
	}
	
	public void initParam() {
		//Mark border type
		HashMap<NodeType, Function> mapNTF = new HashMap<NodeType, Function>();
		mapNTF.put(NodeType.Robin, null);
		mesh.markBorderNode(mapNTF);
		
		//Use element library to assign degree of freedom (DOF) to element
		ElementList eList = mesh.getElementList();
		FELinearTriangle fe = new FELinearTriangle();
		//FEBilinearRectangle fe = new FEBilinearRectangle();
		for(int i=1;i<=eList.size();i++)
			fe.assign(eList.at(i));
	}
	
	public static void plotVector(Mesh mesh, Vector v, String fileName) {
	    MeshWriter writer = new MeshWriter(mesh);
	    if(!outputFolder.isEmpty()) {
		    File file = new File("./"+outputFolder);
			if(!file.exists()) {
				file.mkdirs();
			}
	    }
	    writer.writeTechplot("./"+outputFolder+"/"+fileName, v);
	}
		
	public static void plotFunction(Mesh mesh, Function fun, String fileName) {
	    NodeList list = mesh.getNodeList();
	    int nNode = list.size();
		Variable var = new Variable();
		Vector v = new SparseVector(nNode);
	    for(int i=1;i<=nNode;i++) {
	    	Node node = list.at(i);
	    	var.setIndex(node.globalIndex);
	    	var.set("x", node.coord(1));
	    	var.set("y", node.coord(2));
	    	v.set(i, fun.value(var));
	    }
	    plotVector(mesh,v,fileName);
	}
	
	public Vector solverOneStep(final int step) {
//		f = new AbstractFunction("x","y") {
//			@Override
//			public double value(Variable v) {
//				double x = v.get("x");
//				double y = v.get("y");
//				if(Math.sqrt((x)*(x)+y*y)<0.05)
//					return 10.0*Math.cos(4*Math.PI*Dt*step);
//				//else if(Math.sqrt((x-5)*(x-5)+y*y)<0.3)
//				//	return 100.0;
//				//else if(Math.sqrt((x-5)*(x-5)+(y-5)*(y-5))<0.3)
//				//	return 100.0;
//				else
//					return 0.0;
//			}
//		};
//		System.out.println(f.value(new Variable("x",0).set("y", 0)));
		
		FC k = new FC(Dt*Dt*c2);
		//Function ff = this.u_n1.mult(FC.c(2.0)).minus(this.u_n).plus(f.mult(FC.c(Dt*Dt)));
		Function ff = this.u_n1.X(FC.c(2.0)).M(this.u_n).X(FC.c(1+13*Dt*step));
		weakForm.setF(ff);
		weakForm.setParam(k, FC.c(1.0).X(FC.c(1+13*Dt*step)), null, k);
		
		//Assemble
		AssemblerScalar assembler = new AssemblerScalar(mesh, weakForm);
		System.out.println("Begin Assemble...");
		assembler.assemble();
		Matrix stiff = assembler.getStiffnessMatrix();
		Vector load = assembler.getLoadVector();
		//Boundary condition
		//assembler.imposeDirichletCondition(new FC(0.0));
		System.out.println("Assemble done!");
		
		//CG
		Solver solver = new Solver();
		SparseVector u = new SparseVector(load.getDim(),1.0);
		System.out.println("begin construct AlgebraMatrix...");
		AlgebraMatrix algStiff = new CompressedRowMatrix((SparseMatrix)stiff,true);
		System.out.println("end construct AlgebraMatrix!");
		FullVector algLoad = new FullVector((SparseVector)load);
		FullVector algU = new FullVector(u);
		solver.CG(algStiff, algLoad, algU);
		double[] data = algU.getData();
		for(int i=0;i<data.length;i++) {
			u.set(i+1, data[i]);
		}
		//System.out.println("u=");
		//for(int i=1;i<=u.getDim();i++)
		//	System.out.println(String.format("%.3f", u.get(i)));	
	    
		plotVector(mesh, u, String.format("u_t%03d.dat",step));
		
		this.u_n = this.u_n1;
		this.u_n1 = u_n;

	    return u;
		
	}
	
	public void run() {
		readMesh();
		initParam();
		
		//Time step size
		Dt = 0.004;
		c2 = 500;
		
		u0 = new AbstractFunction("x","y") {
			
//[-3,3]x[-3,3]			
//			@Override
//			public double value(Variable v) {
//				double x = v.get("x");
//				double y = v.get("y");
//				if(Math.sqrt(x*x+y*y)<0.01)
//					return 10.0;
//				else
//					return 0.0;
//			}
			
//凸型区域[-40,40]x[-30,25]			
			@Override
			public double value(Variable v) {
				double x = v.get("x");
				double y = v.get("y");
				if(Math.sqrt(x*x+y*y)<0.55)
					return 20.0;
				else if(Math.sqrt((x-7)*(x-7)+y*y)<0.5)
					return 20.0;
				else if(Math.sqrt((x+15)*(x+15)+y*y)<0.5)
					return 20.0;
				else if(Math.sqrt(x*x+(y-12)*(y-12))<0.2)
					return 100.0;
				else if(Math.sqrt(x*x+(y+10)*(y+10))<2)
					return 5.0;
				else
					return 0.0;
			}
			
//[-100,100]x[-100,100]
//			@Override
//			public double value(Variable v) {
//				double x = v.get("x");
//				double y = v.get("y");
//				if(Math.sqrt((x-30)*(x-30)+y*y)<2)
//					return 100.0;
//				else if(Math.sqrt((x+30)*(x+30)+y*y)<2)
//					return 100.0;
//				else
//					return 0.0;
//			}
		};
		u_n = u0;
		u_n1 = u0;
		for(int i=1;i<=300;i++) {
			Vector rlt = solverOneStep(i);
			u_n1 = new Vector2Function(rlt);
		}
	}
	
	public static void main(String[] args) {
		Wave wave = new Wave();
		wave.run();
	}
}
