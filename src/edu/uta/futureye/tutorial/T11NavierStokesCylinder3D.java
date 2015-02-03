package edu.uta.futureye.tutorial;

import java.util.HashMap;

import edu.uta.futureye.algebra.SparseBlockMatrix;
import edu.uta.futureye.algebra.SparseBlockVector;
import edu.uta.futureye.algebra.SparseVectorHashMap;
import edu.uta.futureye.algebra.intf.SparseVector;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.algebra.solver.SchurComplementStokesSolver;
import edu.uta.futureye.application.DataReader;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.NodeType;
import edu.uta.futureye.function.AbstractMathFun;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.basic.SpaceVectorFunction;
import edu.uta.futureye.function.basic.Vector2Function;
import edu.uta.futureye.function.intf.MathFun;
import edu.uta.futureye.function.intf.VectorFunction;
import edu.uta.futureye.io.MeshReader;
import edu.uta.futureye.lib.assembler.AssemblerVector;
import edu.uta.futureye.lib.element.FETrilinearV_ConstantP;
import edu.uta.futureye.lib.element.FiniteElementType;
import edu.uta.futureye.lib.weakform.WeakFormNavierStokes3D;
import edu.uta.futureye.util.Constant;
import edu.uta.futureye.util.container.ElementList;


/**
 * Problem: 3D Navier-Stokes, flow around a cylinder
 * 
 * Ref.
 * 1. Alexander N. BROOKS Streamline Upwind/Petrow-Galerkin Formulations for Convection Dominated
 *    Flows With Particular Emphasis on the Incompressible Navier-Stokes Equations
 *    
 * 2. M. Schafer and S. Turek Benchmark Computations of Laminar Flow Arond a Cylinder
 * @author liuyueming
 *
 */
public class T11NavierStokesCylinder3D {
	protected String outputFolder = "tutorial/NavierStokesCylinder3D";
	protected String file = null;
	
	protected Mesh mesh = null;
	protected Mesh meshOld = null;
	
	//Navier-Stokes Weak Form (For Picard Iteration)
	protected WeakFormNavierStokes3D weakForm = new WeakFormNavierStokes3D();
	//Assembler
	protected AssemblerVector assembler = null;
	//Dirichlet boundary condition
	protected VectorFunction diri = null;
	//Previous Velocity
	protected VectorFunction U = new SpaceVectorFunction(3);
	
	//delta t
	protected double dt = 0.02;
	//viscosity
	protected double nu = 0.001; //0.0005
	
	FiniteElementType fe = null;
	
	int maxNonlinearIter = 30;
	double nonlinearError = 1e-2;
	int maxTimeStep = 1000;
	
	
	/**
	 * 
	 * @param testCaseNo
	 */
	public void init(int testCaseNo) {
		//Read mesh from input file
		//file = "benchmark_cylinder3_3D"; //rectangle mesh
		file = "benchmark_cylinder31_3D"; //rectangle mesh
		
		MeshReader reader = new MeshReader(file+".grd");
		MeshReader reader2 = new MeshReader(file+".grd");
		mesh = reader.read3DMesh();
		meshOld = reader2.read3DMesh();
		mesh.nVertex = mesh.getNodeList().size();
		
		//Geometry relationship
		mesh.computeNodeBelongsToElements();
		mesh.computeGlobalEdge();
		
		ElementList eList = mesh.getElementList();
//		NodeList nodes = mesh.getNodeList();
//		for(int i=1;i<=eList.size();i++) {
//			System.out.println(i+"  " + eList.at(i));
//		}
		
		Vector v = new SparseVectorHashMap(mesh.getNodeList().size());
		Tools.plotVector(mesh, outputFolder, "test.dat", v);

		//Mark boundary type of u
		HashMap<NodeType, MathFun> mapNTF_u = new HashMap<NodeType, MathFun>();
		//Dirichlet boundary of u
		mapNTF_u.put(NodeType.Dirichlet, new AbstractMathFun("x","y","z") {
			@Override
			public double apply(Variable v) {
				double x = v.get("x");
				double y = v.get("y");
				double z = v.get("z");
				double redius = 0.05;
				double H = 0.41;
				double shift_x = -0.3;
				//upside and down side boundary
				if(y < Constant.meshEps || Math.abs(y-H) < Constant.meshEps)
					return 1;
				//back and front side boundary
				if(z < Constant.meshEps || Math.abs(z-H) < Constant.meshEps)
					return 1;
				//left side boundary
				if(Math.abs(x-shift_x) < Constant.meshEps)
					return 1;
				//cylinder boundary
				if(Math.sqrt((x-0.2)*(x-0.2)+(y-0.2)*(y-0.2)) <= redius+Constant.meshEps)
					return 1;
				return 0;
			}
		});
		//Neumann boundary of u
		mapNTF_u.put(NodeType.Neumann, null);
		
		//Mark boundary type of p
		HashMap<NodeType, MathFun> mapNTF_p = new HashMap<NodeType, MathFun>();
		//Dirichlet boundary of p
		mapNTF_p.put(NodeType.Dirichlet, new AbstractMathFun("x","y","z") {
			@Override
			public double apply(Variable v) {
				double x = v.get("x");
				//right side p=0
				if(Math.abs(x-2.2) < Constant.meshEps)
					return 1;
				return 0;
			}
		});
		//Neumann boundary of p
		mapNTF_p.put(NodeType.Neumann, null);
		
		mesh.markBorderNode(1,mapNTF_u); //u's boundary type
		mesh.markBorderNode(2,mapNTF_u); //v's boundary type same as u
		mesh.markBorderNode(3,mapNTF_u); //w's boundary type same as u
		mesh.markBorderNode(4,mapNTF_p);
		
		mesh.writeNodesInfo(String.format("./%s/%s__MeshInfo_u.dat",outputFolder,file), 1);
		mesh.writeNodesInfo(String.format("./%s/%s__MeshInfo_v.dat",outputFolder,file), 2);
		mesh.writeNodesInfo(String.format("./%s/%s__MeshInfo_w.dat",outputFolder,file), 3);
		mesh.writeNodesInfo(String.format("./%s/%s__MeshInfo_p.dat",outputFolder,file), 4);

		//Use element library to assign degree of freedom (DOF) to element
		fe = new FETrilinearV_ConstantP();
		fe.initDOFIndexGenerator(mesh);
		for(int i=1;i<=eList.size();i++) {
			fe.assignTo(eList.at(i));
//			eList.at(i).printDOFInfo();
		}


	}
	
	/**
	 * u_max=0.45 at left side boundary, other u=v=w=0
	 * p=0 at right side boundary
	 * 
	 * Result Re=20
	 */
	public void imposeDirichletCondition1() {
		diri = new SpaceVectorFunction(4);
		diri.set(1, new AbstractMathFun("x","y","z") {
					@Override
					public double apply(Variable v) {
						double x = v.get("x");
						double y = v.get("y");
						double z = v.get("z");
						double H  = 0.41;
						double Um = 0.45; //!!!
						double shift_x = -0.3;
						//left side boundary
						if(Math.abs(x-shift_x) < Constant.meshEps)
							return 16*Um*y*z*(H-y)*(H-z)/(H*H*H*H);
						else
							return 0.0;
					}
				});
		diri.set(2, FC.C0);
		diri.set(3, FC.C0);
		diri.set(4, FC.C0);		
	}
	
	/**
	 * u_max=2.25 at left side boundary, other u=v=w=0
	 * p=0 at right side boundary
	 * 
	 * Result Re=100
	 */
	public void imposeDirichletCondition2() {
		diri = new SpaceVectorFunction(4);
		diri.set(1, new AbstractMathFun("x","y","z") {
			@Override
			public double apply(Variable v) {
				double x = v.get("x");
				double y = v.get("y");
				double z = v.get("z");
				double H  = 0.41;
				double Um = 2.25; //!!!
				double shift_x = -0.3;
				//left side boundary
				if(Math.abs(x-shift_x) < Constant.meshEps)
					return 16*Um*y*z*(H-y)*(H-z)/(H*H*H*H);
				else
					return 0.0;
					}
				});
		diri.set(2, FC.C0);
		diri.set(3, FC.C0);
		diri.set(4, FC.C0);		
	}

	/**
	 * u_max=2.25 at left side boundary, u=1 at upside,down,back and front side, other u=v=w=0
	 * p=0 at right side boundary
	 * 
	 * Result Re=100
	 */
	public void imposeDirichletCondition3() {
		diri = new SpaceVectorFunction(4);
		diri.set(1, new AbstractMathFun("x","y","z") {
					@Override
					public double apply(Variable v) {
						double x = v.get("x");
						double y = v.get("y");
						double z = v.get("z");
						double H  = 0.41;
						double Um = 2.25;//0.45; //!!!
						double shift_x = -0.3;
						//upside and down side boundary
						if(y < Constant.meshEps || Math.abs(y-0.41) < Constant.meshEps)
							return 0.3;
						//back and front side boundary
						if(z < Constant.meshEps || Math.abs(z-0.41) < Constant.meshEps)
							return 0.3;
						//left side boundary
						if(Math.abs(x-shift_x) < Constant.meshEps)
							return 16*Um*(y-0.03)*z*(H-(y-0.03))*(H-z)/(H*H*H*H);
						else
							return 0.0;
					}
				});
		diri.set(2, FC.C0);
		diri.set(3, FC.C0);
		diri.set(4, FC.C0);	
	}
	
	public SparseBlockVector nonlinearIter(int time, int nIter, SpaceVectorFunction uk) {
		//Right hand side(RHS): f = (0,0)'
		if(time==0)
			weakForm.setF(new SpaceVectorFunction(FC.C0,FC.C0,FC.C0));
		else
			weakForm.setF(new SpaceVectorFunction(
					uk.get(1).D(dt),
					uk.get(2).D(dt),
					uk.get(3).D(dt)
					));
			
		weakForm.setParam(FC.c(nu),U,FC.c(1.0/dt));
		
		assembler = new AssemblerVector(mesh, weakForm,fe);
		assembler.assemble();
		SparseBlockMatrix stiff = assembler.getStiffnessMatrix();
		SparseBlockVector load = assembler.getLoadVector();

		assembler.imposeDirichletCondition(diri);
		
		SchurComplementStokesSolver solver = 
			new SchurComplementStokesSolver(stiff,load);
		//solver.setCGInit(0.5);
		//solver.debug = true;
		return solver.solve3D();
		
	}
	
	public SparseBlockVector nonlinearIterSteady(int nIter, SpaceVectorFunction uk) {
		weakForm.setF(new SpaceVectorFunction(FC.C0,FC.C0,FC.C0));
		weakForm.setParam(FC.c(nu),U,FC.C0);
		
		assembler = new AssemblerVector(mesh, weakForm,fe);
		System.out.println("Begin Assemble...");
		assembler.assemble();
		SparseBlockMatrix stiff = assembler.getStiffnessMatrix();
		SparseBlockVector load = assembler.getLoadVector();
		System.out.print("imposeDirichletCondition...");
		long begin = System.currentTimeMillis();
		assembler.imposeDirichletCondition(diri);
		long end = System.currentTimeMillis();
		System.out.println((end-begin)+"ms");
		System.out.println("Assemble done!");
		
		SchurComplementStokesSolver solver = 
			new SchurComplementStokesSolver(stiff,load);
		//solver.setCGInit(0.5);
		//solver.debug = true;
		return solver.solve3D();
		
	}
	
	public void run(int startTimeStep, int testCaseNo, boolean bSteady) {
		init(testCaseNo);
		
		if(!bSteady && startTimeStep>0) {
			Vector vecU = DataReader.readVector(String.format("./%s/%s_uvw_final_t%02d.dat",
					outputFolder,file,startTimeStep),3);
			Vector vecV = DataReader.readVector(String.format("./%s/%s_uvw_final_t%02d.dat",
					outputFolder,file,startTimeStep),4);
			Vector vecW = DataReader.readVector(String.format("./%s/%s_uvw_final_t%02d.dat",
					outputFolder,file,startTimeStep),5);
			U.set(1, new Vector2Function(vecU));
			U.set(2, new Vector2Function(vecV));
			U.set(3, new Vector2Function(vecW));
		} else if(bSteady && startTimeStep>0){ //startTimeStep will be iteration restart step in steady case
			Vector vecU = DataReader.readVector(String.format("./%s/%s_uvw_steady_%02d.dat",
					outputFolder,file,startTimeStep),3);
			Vector vecV = DataReader.readVector(String.format("./%s/%s_uvw_steady_%02d.dat",
					outputFolder,file,startTimeStep),4);
			Vector vecW = DataReader.readVector(String.format("./%s/%s_uvw_steady_%02d.dat",
					outputFolder,file,startTimeStep),5);
			U.set(1, new Vector2Function(vecU));
			U.set(2, new Vector2Function(vecV));
			U.set(3, new Vector2Function(vecW));
		} else {
			U.set(1, FC.C0);
			U.set(2, FC.C0);
			U.set(3, FC.C0);
		}
		SparseBlockVector u = null;
		if(bSteady) System.out.println(">>>>>>>>>>>>>>>>>>>steady");
		SpaceVectorFunction uk = new SpaceVectorFunction(3);
		for(int time=startTimeStep+1;time<this.maxTimeStep;time++) {
			if(!bSteady) System.out.println(">>>>>>>>>>>>>>>>>>>time="+time);
			uk.set(1, U.get(1));
			uk.set(2, U.get(2));
			uk.set(3, U.get(3));
			for(int iter=0;iter<this.maxNonlinearIter;iter++) {
				if(bSteady)
					u = nonlinearIterSteady(iter, uk);
				else
					u = nonlinearIter(time, iter, uk);
				
				//Compute norm of delta_u (not delta_v)
				int dim = u.getBlock(1).getDim();
				SparseVector delta_u = new SparseVectorHashMap(dim);
				for(int i=1;i<=dim;i++)
					delta_u.set(i, 
							u.getBlock(1).get(i)-
							U.get(1).apply(new Variable().setIndex(i)));
				
				U.set(1, new Vector2Function(u.getBlock(1)));
				U.set(2, new Vector2Function(u.getBlock(2)));
				U.set(3, new Vector2Function(u.getBlock(3)));
	
				System.out.println("Iter="+iter+" Error Norm2 (||u1_k+1 - u1_k||) = "+delta_u.norm2());
				
				if(delta_u.norm2() < this.nonlinearError) {
					String s = "_t%03d";
					if(bSteady) s = "_steady";
					Tools.plotVector(mesh, outputFolder, String.format("%s_uvw_final"+s+".dat",file,time), 
							u.getBlock(1), u.getBlock(2), u.getBlock(3));
					Tools.plotVector(meshOld, outputFolder, String.format("%s_p_final"+s+".dat",file,time), 
							Tools.valueOnElement2Node(mesh,u.getBlock(4)));
					if(bSteady)
						return;
					else
						break;
				} else {
					if(bSteady) {
						Tools.plotVector(mesh, outputFolder, String.format("%s_uvw_steady_%02d.dat",file,time+iter), 
								u.getBlock(1), u.getBlock(2), u.getBlock(3));
						Tools.plotVector(meshOld, outputFolder, String.format("%s_p_steady_%02d.dat",file,time+iter), 
								Tools.valueOnElement2Node(mesh,u.getBlock(4)));
					} else {
						Tools.plotVector(mesh, outputFolder, String.format("%s_uvw%02d_%02d.dat",file,time,iter), 
								u.getBlock(1), u.getBlock(2), u.getBlock(3));
						Tools.plotVector(meshOld, outputFolder, String.format("%s_p%02d_%02d.dat",file,time,iter), 
								Tools.valueOnElement2Node(mesh,u.getBlock(4)));
					}
				}
			}
			if(bSteady) return;
		}
	}

	
	/**
	 * args[0]: test case number
	 * args[1]: is steady
	 * args[2]: delta t
	 * args[3]: start time step
	 * args[4]: max time step

	 * @param args
	 */
	public static void main(String[] args) {
		T11NavierStokesCylinder3D NS = new T11NavierStokesCylinder3D();
		int testCaseNo = 3;
		boolean bSteady = false;
		NS.dt = 0.02;
		int startTimeStep = 0;
		if(args.length >= 1)
			testCaseNo = Integer.parseInt(args[0]);
		if(args.length >= 2)
			bSteady = Boolean.parseBoolean(args[1]);
		if(args.length >= 3)
			NS.dt = Double.parseDouble(args[2]);
		if(args.length >= 4)
			startTimeStep = Integer.parseInt(args[3]);
		if(args.length >= 5)
			NS.maxTimeStep = Integer.parseInt(args[4]);

		System.out.println("testCaseNo="+testCaseNo);
		System.out.println("bSteady="+bSteady);
		System.out.println("dt="+NS.dt);
		System.out.println("startTimeStep="+startTimeStep);
		System.out.println("maxTimeStep="+NS.maxTimeStep);
		System.out.println("------default values------");
		System.out.println("mu="+NS.nu);
		System.out.println("maxNonlinearIter="+NS.maxNonlinearIter);
		System.out.println("nonlinearError="+NS.nonlinearError);
		
		//NS.imposeDirichletCondition1();
		//NS.imposeDirichletCondition2();
		NS.imposeDirichletCondition3();
		NS.run(startTimeStep,testCaseNo,bSteady);
	}
}