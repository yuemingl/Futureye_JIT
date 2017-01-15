package edu.uta.futureye.tutorial;

import java.util.HashMap;

import edu.uta.futureye.algebra.SparseBlockMatrix;
import edu.uta.futureye.algebra.SparseBlockVector;
import edu.uta.futureye.algebra.SparseVectorHashMap;
import edu.uta.futureye.algebra.intf.SparseVector;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.algebra.solver.SchurComplementStokesSolver;
import edu.uta.futureye.application.DataReader;
import edu.uta.futureye.core.EdgeLocal;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.core.NodeLocal;
import edu.uta.futureye.core.NodeType;
import edu.uta.futureye.core.Vertex;
import edu.uta.futureye.function.MultiVarFunc;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.basic.SpaceVectorFunction;
import edu.uta.futureye.function.basic.Vector2Function;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.function.intf.VectorMathFunc;
import edu.uta.futureye.io.MeshReader;
import edu.uta.futureye.lib.assembler.AssemblerVector;
import edu.uta.futureye.lib.element.FEBilinearV_ConstantP;
import edu.uta.futureye.lib.element.FEQuadraticV_LinearP;
import edu.uta.futureye.lib.element.FiniteElementType;
import edu.uta.futureye.lib.weakform.WeakFormNavierStokes2D;
import edu.uta.futureye.util.Constant;
import edu.uta.futureye.util.FutureyeException;
import edu.uta.futureye.util.container.ElementList;
import edu.uta.futureye.util.container.ObjIndex;
import edu.uta.futureye.util.container.ObjList;
import static edu.uta.futureye.function.FMath.*;


/**
 * Problem: 2D Navier-Stokes, lid-driven cavity
 *
 * Ref.
 * 1. T.E. Teaduyar Stabilized Finite Element Formulations for Incompressible Flow Computations
 * 
 * 2. Rajiv Sampath and Nicholas Zabaras Design and object-oriented implementation of a 
 *    preconditioned-stabilized incompressible NaiverStokes solver using equal-order-interpolation 
 *    velocity pressure elements
 *
 * @author liuyueming
 *
 */
public class T11NavierStokesBox {
	protected String outputFolder = "tutorial/NavierStokesBox";
	protected String file = null;
	
	protected Mesh mesh = null;
	protected Mesh meshOld = null;
	
	//Navier-Stokes Weak Form (For Picard Iteration)
	protected WeakFormNavierStokes2D weakForm = new WeakFormNavierStokes2D();
	//Assembler
	protected AssemblerVector assembler = null;
	//Dirichlet boundary condition
	protected VectorMathFunc diri = null;
	//Previous Velocity
	protected VectorMathFunc U = new SpaceVectorFunction(2);
	
	//delta t
	protected double dt = 0.02;
	//viscosity
	protected double nu = 0.001; 
	
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
		if(testCaseNo == 1)
			file = "stokes_box";//[-3,3]*[-3,3]
		else if(testCaseNo == 2) 
			file = "stokes_box2";//[0,1]*[0,1] 64*64
		else if(testCaseNo == 3) 
			file = "stokes_box3";//[0,1]*[0,1] 32*32
		else
			throw new FutureyeException("testCaseNo should be 1,2");
		
		MeshReader reader = new MeshReader(file+".grd");
		MeshReader reader2 = new MeshReader(file+".grd");
		mesh = reader.read2DMesh();
		meshOld = reader2.read2DMesh();
		mesh.nVertex = mesh.getNodeList().size();
		
		//Add nodes for quadratic element, original mesh stored in meshOld
		if(testCaseNo == 1) {
			for(int i=1;i<=mesh.getElementList().size();i++) {
				Element e = mesh.getElementList().at(i);
				e.adjustVerticeToCounterClockwise();
				ObjList<EdgeLocal> edges = e.edges();
				int nNode = e.nodes.size();
				for(int j=1;j<=edges.size();j++) {
					EdgeLocal edge = edges.at(j);
					Vertex l = edge.beginVertex();
					Vertex r = edge.endVertex();
					double cx = (l.coord(1)+r.coord(1))/2.0;
					double cy = (l.coord(2)+r.coord(2))/2.0;
					Node node = new Node(mesh.getNodeList().size()+1, cx,cy);
					Node findNode = mesh.findNode(node);
					if(findNode == null) {
						edge.addEdgeNode(new NodeLocal(++nNode,node));
						mesh.addNode(node);
					} else {
						edge.addEdgeNode(new NodeLocal(++nNode,findNode));
					}
				}
				e.applyChange();
			}
		}
		//Geometry relationship
		mesh.computeNodeBelongsToElements();
		
		ElementList eList = mesh.getElementList();
//		NodeList nodes = mesh.getNodeList();
//		for(int i=1;i<=eList.size();i++) {
//			System.out.println(i+"  " + eList.at(i));
//		}
		
		//Mark border type
		HashMap<NodeType, MathFunc> mapNTF_uv = new HashMap<NodeType, MathFunc>();
		mapNTF_uv.put(NodeType.Dirichlet, null);
		
		HashMap<NodeType, MathFunc> mapNTF_p = new HashMap<NodeType, MathFunc>();
		mapNTF_p.put(NodeType.Neumann, null);
		
		mesh.markBorderNode(new ObjIndex(1,2),mapNTF_uv);
		mesh.markBorderNode(3,mapNTF_p);

		//Use element library to assign degree of freedom (DOF) to element
		if(testCaseNo == 1)
			fe = new FEQuadraticV_LinearP();
		else if(testCaseNo == 2 || testCaseNo == 3)
			fe = new FEBilinearV_ConstantP();
		fe.initDOFIndexGenerator(mesh);
		for(int i=1;i<=eList.size();i++) {
			fe.assignTo(eList.at(i));
			//eList.at(i).printDOFInfo();
		}

		diri = new SpaceVectorFunction(3);
		diri.set(1, new MultiVarFunc("", "x", "y") {
//					@Override
//					public double apply(Variable v) {
//						double y = v.get("y");
//						if(Math.abs(y-1.0)<Constant.meshEps)
//							return 1.0;
//						else
//							return 0.0;
//					}

					@Override
					public double apply(double... args) {
						double y = args[this.argIdx[1]];
						if(Math.abs(y-1.0)<Constant.meshEps)
							return 1.0;
						else
							return 0.0;
					}
				});
		diri.set(2, C0);
		diri.set(3, C0);
		
	}
	
	public SparseBlockVector nonlinearIter(int time, int nIter, SpaceVectorFunction uk) {
		//Right hand side(RHS): f = (0,0)'
		if(time==0)
			weakForm.setF(new SpaceVectorFunction(C0,C0));
		else
			weakForm.setF(new SpaceVectorFunction(uk.get(1).D(dt),uk.get(2).D(dt)));
		
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
		return solver.solve2D();
		
	}	
	
	public SparseBlockVector nonlinearIterSteady(int nIter, SpaceVectorFunction uk) {
		//Right hand side(RHS): f = (0,0)'
		weakForm.setF(new SpaceVectorFunction(C0,C0));
		weakForm.setParam(FC.c(nu),U,C0);
		
		assembler = new AssemblerVector(mesh, weakForm,fe);
		assembler.assemble();
		SparseBlockMatrix stiff = assembler.getStiffnessMatrix();
		SparseBlockVector load = assembler.getLoadVector();
		System.out.print("imposeDirichletCondition...");
		long begin = System.currentTimeMillis();
		assembler.imposeDirichletCondition(diri);
		long end = System.currentTimeMillis();
		System.out.println((end-begin)+"ms");
		
		SchurComplementStokesSolver solver = 
			new SchurComplementStokesSolver(stiff,load);
		//solver.setCGInit(0.5);
		//solver.debug = true;
		return solver.solve2D();
		
	}	
	
	public void run(int startTimeStep, int testCaseNo, boolean bSteady) {
		init(testCaseNo);
		if(bSteady) startTimeStep=0;
		
		if(startTimeStep>0) {
			Vector vecU = DataReader.readVector(String.format("./%s/%s_uv_final_t%03d.dat",
					outputFolder,file,startTimeStep),3);
			Vector vecV = DataReader.readVector(String.format("./%s/%s_uv_final_t%03d.dat",
					outputFolder,file,startTimeStep),4);
			U.set(1, new Vector2Function(vecU));
			U.set(2, new Vector2Function(vecV));
		} else {
			U.set(1, C0);
			U.set(2, C0);
		}
		
		
		SparseBlockVector u = null;
		if(bSteady) System.out.println(">>>>>>>>>>>>>>>>>>>steady");
		SpaceVectorFunction uk = new SpaceVectorFunction(2);
		for(int time=startTimeStep+1;time<this.maxTimeStep;time++) {
			if(!bSteady) System.out.println(">>>>>>>>>>>>>>>>>>>time="+time);

			uk.set(1, U.get(1));
			uk.set(2, U.get(2));
			for(int iter=0;iter<this.maxNonlinearIter;iter++) {
				if(bSteady)
					u = nonlinearIterSteady(iter, uk);
				else
					u = nonlinearIter(time, iter, uk);
				
				//Compute norm of delta_u (not including delta_v)
				int dim = u.getBlock(1).getDim();
				SparseVector delta_u = new SparseVectorHashMap(dim);
				for(int i=1;i<=dim;i++)
					delta_u.set(i, 
							
							u.getBlock(1).get(i)-
							U.get(1).apply(new Variable().setIndex(i)));
				
				U.set(1, new Vector2Function(u.getBlock(1)));
				U.set(2, new Vector2Function(u.getBlock(2)));
	
				System.out.println("Iter="+iter+" Error Norm2 (||u1_k+1 - u1_k||) = "+delta_u.norm2());
				
				if(delta_u.norm2() < this.nonlinearError) {
					String s = "_t%03d";
					if(bSteady) s = "_steady";
					Tools.plotVector(mesh, outputFolder, String.format("%s_uv_final"+s+".dat",file,time), 
							u.getBlock(1), u.getBlock(2));
					Tools.plotVector(meshOld, outputFolder, String.format("%s_p_final"+s+".dat",file,time), 
							Tools.valueOnElement2Node(mesh,u.getBlock(3)));
					if(bSteady)
						return;
					else
						break;
				} else {
					Tools.plotVector(mesh, outputFolder, String.format("%s_uv%02d_%02d.dat",file,time,iter), 
							u.getBlock(1), u.getBlock(2));
					Tools.plotVector(meshOld, outputFolder, String.format("%s_p%02d_%02d.dat",file,time,iter), 
							Tools.valueOnElement2Node(mesh, u.getBlock(3)));
				}
			}
		}
	}
	
	/**
	 * args[0]: test case number
	 * args[1]: is steady
	 * args[2]: delta t
	 * args[3]: start time step
	 * args[4]: max start time step
	 * 
	 * @param args
	 */
	public static void main(String[] args) {
		T11NavierStokesBox NSB = new T11NavierStokesBox();
		NSB.dt = 0.02;
		int startTimeStep = 0;
		int testCaseNo = 3;
		boolean bSteady = false;
		if(args.length >= 1)
			testCaseNo = Integer.parseInt(args[0]);
		if(args.length >= 2)
			bSteady = Boolean.parseBoolean(args[1]);
		if(args.length >= 3)
			NSB.dt = Double.parseDouble(args[2]);
		if(args.length >= 4)
			startTimeStep = Integer.parseInt(args[3]);
		if(args.length >= 5)
			NSB.maxTimeStep = Integer.parseInt(args[4]);
		
		System.out.println("testCaseNo="+testCaseNo);
		System.out.println("bSteady="+bSteady);
		System.out.println("dt="+NSB.dt);
		System.out.println("timeStep="+startTimeStep);
		System.out.println("maxTimeStep="+NSB.maxTimeStep);
		System.out.println("default values:");
		System.out.println("mu="+NSB.nu);
		System.out.println("maxNonlinearIter="+NSB.maxNonlinearIter);
		System.out.println("nonlinearError="+NSB.nonlinearError);
		
		NSB.run(startTimeStep,testCaseNo,bSteady);
	}
}
