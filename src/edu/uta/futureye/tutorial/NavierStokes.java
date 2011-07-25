package edu.uta.futureye.tutorial;

import java.util.HashMap;

import edu.uta.futureye.algebra.SchurComplementStokesSolver;
import edu.uta.futureye.algebra.SparseVector;
import edu.uta.futureye.algebra.intf.BlockMatrix;
import edu.uta.futureye.algebra.intf.BlockVector;
import edu.uta.futureye.core.EdgeLocal;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.core.NodeLocal;
import edu.uta.futureye.core.NodeType;
import edu.uta.futureye.core.Vertex;
import edu.uta.futureye.function.AbstractFunction;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.basic.SpaceVectorFunction;
import edu.uta.futureye.function.basic.Vector2Function;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.function.intf.VectorFunction;
import edu.uta.futureye.io.MeshReader;
import edu.uta.futureye.lib.assembler.AssemblerVector;
import edu.uta.futureye.lib.element.FEQuadraticV_LinearP;
import edu.uta.futureye.lib.weakform.WeakFormNavierStokes;
import edu.uta.futureye.util.Constant;
import edu.uta.futureye.util.container.ElementList;
import edu.uta.futureye.util.container.NodeList;
import edu.uta.futureye.util.container.ObjIndex;
import edu.uta.futureye.util.container.ObjList;


/**
 * Problem: Navier-Stokes
 * 
 * @author liuyueming
 *
 */
public class NavierStokes {
	protected String file = "benchmark_cylinder1";
	protected static String outputFolder = "tutorial\\NavierStokes";
	protected Mesh mesh = null;
	protected Mesh meshOld = null;
	//Quadratic Velocity - Linear Pressure Element
	protected FEQuadraticV_LinearP fe = new FEQuadraticV_LinearP();
	//Stokes weak form
	protected WeakFormNavierStokes weakForm = new WeakFormNavierStokes();
	//Assembler
	protected AssemblerVector assembler = null;
	//Dirichlet boundary condition
	protected VectorFunction diri = null;
	//Previous velocity
	protected VectorFunction U = new SpaceVectorFunction(2);
	
	public void init() {
		//Read a triangle mesh from an input file
		MeshReader reader = new MeshReader(file+".grd");
		MeshReader reader2 = new MeshReader(file+".grd");
		mesh = reader.read2DMesh();
		meshOld = reader2.read2DMesh();
		mesh.nVertex = mesh.getNodeList().size();
		
		//Add nodes for quadratic element
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
				Node findNode = mesh.containNode(node);
				if(findNode == null) {
					edge.addEdgeNode(new NodeLocal(++nNode,node));
					mesh.addNode(node);
				} else {
					edge.addEdgeNode(new NodeLocal(++nNode,findNode));
				}
			}
			e.applyChange();
		}
		
		//Geometry relationship
		mesh.computeNodeBelongsToElements();
		
		ElementList eList = mesh.getElementList();
		NodeList nodes = mesh.getNodeList();
//		for(int i=1;i<=eList.size();i++) {
//			System.out.println(i+"  " + eList.at(i));
//		}
		
		//Mark border type
		//指定u,v边界
		HashMap<NodeType, Function> mapNTF_uv = new HashMap<NodeType, Function>();
		mapNTF_uv.put(NodeType.Dirichlet, new AbstractFunction("x","y") {
			@Override
			public double value(Variable v) {
				double x = v.get("x");
				double y = v.get("y");
				//upside and down side
				if(y < Constant.meshEps || Math.abs(y-0.41) < Constant.meshEps)
					return 1;
				//left side
				else if(x < Constant.meshEps)
					return 1;
				
				//cylinder
				if(Math.sqrt((x-0.2)*(x-0.2)+(y-0.2)*(y-0.2)) <= 0.05+Constant.meshEps)
					return 1;
				
				return 0;
			}
		});
		//u,v其他边界
		mapNTF_uv.put(NodeType.Neumann, null);
		
		//指定p边界
		HashMap<NodeType, Function> mapNTF_p = new HashMap<NodeType, Function>();
		mapNTF_p.put(NodeType.Dirichlet, new AbstractFunction("x","y") {
			@Override
			public double value(Variable v) {
				double x = v.get("x");
				if(Math.abs(x-2.2) < Constant.meshEps)
					return 1;
				else
					return 0;
			}
		});
		//p其他边界
		mapNTF_p.put(NodeType.Neumann, null);
		
		mesh.markBorderNode(new ObjIndex(1,2),mapNTF_uv);
		mesh.markBorderNode(3,mapNTF_p);

		//Use element library to assign degree of freedom (DOF) to element
		for(int i=1;i<=eList.size();i++) {
			System.out.println(i+"  " + eList.at(i));
		}

		fe.initDOFIndexGenerator(nodes.size());
		for(int i=1;i<=eList.size();i++) {
			fe.assignTo(eList.at(i));
			//eList.at(i).printDOFInfo();
		}

		//Boundary condition
		diri = new SpaceVectorFunction(3);
		diri.set(1, new AbstractFunction("x","y") {
					@Override
					public double value(Variable v) {
						double x = v.get("x");
						double y = v.get("y");
						
						double H  = 0.41;
						double Um = 0.3;
						if(x < Constant.meshEps)
							return 4*Um*y*(H-y)/(H*H);
						else
							return 0.0;
					}
				});
		diri.set(2, FC.c0);
		diri.set(3, FC.c0);
	}
	
	public BlockVector nonlinearIter(int nIter) {
		//Right hand side(RHS): f = (0,0)'
		weakForm.setF(new SpaceVectorFunction(FC.c0,FC.c0));
		
		weakForm.setParam(FC.c(0.1),U,FC.c0);
		
		//Robin:  k*u_n + d*u - p\mathbf{n} = 0
		VectorFunction d = new SpaceVectorFunction(2);
		d.set(1, FC.c0);
		d.set(2, FC.c0);
		weakForm.setRobin(d);
		
		assembler = new AssemblerVector(mesh, weakForm,fe);
		System.out.println("Begin Assemble...");
		assembler.assemble();
		BlockMatrix stiff = (BlockMatrix)assembler.getStiffnessMatrix();
		BlockVector load = (BlockVector)assembler.getLoadVector();
		//stiff.print();
		//load.print();
		//System.out.println(load.norm2());

		assembler.imposeDirichletCondition(diri);
		//load.getBlock(1).print();
		//load.getBlock(2).print();
		//load.getBlock(3).print();
		System.out.println("Assemble done!");
		
		SchurComplementStokesSolver solver = 
			new SchurComplementStokesSolver(stiff,load);
		
		return solver.solve();
		
	}
	
	public void run() {
		init();
		
		U.set(1, FC.c0);
		U.set(2, FC.c0);
		BlockVector u = null;
		for(int iter=0;iter<15;iter++) {
			u = nonlinearIter(iter);
			
			int dim = u.getBlock(1).getDim();
			SparseVector tmp = new SparseVector(dim);
			for(int i=1;i<=dim;i++)
				tmp.set(i, 
						u.getBlock(1).get(i)-
						U.get(1).value(new Variable().setIndex(i)));
			
			U.set(1, new Vector2Function(u.getBlock(1)));
			U.set(2, new Vector2Function(u.getBlock(2)));

			System.out.println("u=");
			for(int i=1;i<=u.getDim();i++)
				System.out.println(String.format("%.3f", u.get(i)));	
			Tools.plotVector(mesh, outputFolder, String.format("%s_uv%02d.dat",file,iter), 
					u.getBlock(1), u.getBlock(2));
			Tools.plotVector(meshOld, outputFolder, String.format("%s_p%02d.dat",file,iter), 
					u.getBlock(3));
			
			System.out.println("Error Norm = "+tmp.norm2());
		
		}
	}
	
	public static void main(String[] args) {
		NavierStokes NS = new NavierStokes();
		NS.run();
	}
}
