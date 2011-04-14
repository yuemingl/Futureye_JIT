package edu.uta.futureye.tutorial;

import java.util.HashMap;

import edu.uta.futureye.algebra.SchurComplementStokesSolver;
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
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.function.intf.VectorFunction;
import edu.uta.futureye.io.MeshReader;
import edu.uta.futureye.lib.assembler.AssemblerVector;
import edu.uta.futureye.lib.element.FEQuadraticV_LinearP;
import edu.uta.futureye.lib.weakform.WeakFormStokes;
import edu.uta.futureye.util.Constant;
import edu.uta.futureye.util.container.ElementList;
import edu.uta.futureye.util.container.NodeList;
import edu.uta.futureye.util.container.ObjIndex;
import edu.uta.futureye.util.container.ObjList;


/**
 * Problem:
 *   -\Nabla{k*\Nabla{\mathbf{u}} + \Nabla{p} = \mathbf{f}
 *   div{\mathbf{u}} = 0
 * 
 * Each dim:
 *   -k*(u1_xx+u1_yy) + p_x = f1
 *   -k*(u2_xx+u2_yy) + p_y = f2
 *   u1_x+u2_y              = 0

 * Weak form:
 *   find \mathbf{u} \in H_0^1(div;\Omega), p \in L_2(\Omega)
 *   such that, for all \mathbf{v} \in H_0^1(div;\Omega), q \in L_2(\Omega)
 *   
 *   (\Nabla{\mathbf{v}},k*\Nabla{\mathbf{u}}) - (div{\mathbf{v}},p) 
 *                   + (q,div{\mathbf{u}}) = (\mathbf{v},\mathbf{f})
 *
 *   (v1_x,k*u1_x) + (v1_y,k*u1_y) + (v2_x,k*u2_x) + (v2_y,k*u2_y) 
 *                   - (v1_x+v2_y,p) + (q,u1_x+u2_y) = (v1*f1+v2*f2)      
 *
 * where
 *   \mathbf{u}=(u1,u2): velocity vector field    
 *   \mathbf{f}=(f1,f2): body force
 *   
 * @author liuyueming
 *
 */
public class Stokes {
	public static String outputFolder = "tutorial\\Stokes";
	
	public static void cylinder() {
		//Read a triangle mesh from an input file
//		String file = "cylinder1";
//		String file = "cylinder2";
//		String file = "cylinder2_test";
//		String file = "cylinder4";
//		String file = "cylinder_patch";
		String file = "u_shape";
		
		MeshReader reader = new MeshReader(file+".grd");
		MeshReader reader2 = new MeshReader(file+".grd");
		Mesh mesh = reader.read2DMesh();
		Mesh meshOld = reader2.read2DMesh();
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
		
		for(int i=1;i<=eList.size();i++) {
			System.out.println(i+"  " + eList.at(i));
		}
		
		//Mark border type
		HashMap<NodeType, Function> mapNTF_uv = new HashMap<NodeType, Function>();
		//指定边界
		//mapNTF.put(NodeType.Dirichlet,null);
		mapNTF_uv.put(NodeType.Dirichlet, new AbstractFunction("x","y") {
			@Override
			public double value(Variable v) {
				double x = v.get("x");
				double y = v.get("y");
				//all
				if(Math.abs(y-2.5)<Constant.eps|| Math.abs(y+2.5)<Constant.eps)
					return 1;
				
				//cylinder1 and cylinder2 and cylinder4
//				else if(Math.abs(x+5.0)<0.01)
//						return 1;
				
				//u_shape
				else if(Math.abs(x+5.0)<0.01 && y>0)
						return 1;
				else if(Math.abs(x-5.0)<0.01)
					return 1;
				else if(Math.abs(x-4.0)<0.01)
					return 1;
				else if(Math.abs(y-0.5)<Constant.eps|| Math.abs(y+0.5)<Constant.eps)
					return 1;
				
				//cylinder1
//				else if(Math.sqrt(x*x+y*y)<(0.5+0.1))
//					return 1;
				
				//cylinder2 and cylinder4
//				else if(Math.sqrt(x*x+(y-1)*(y-1))<(0.5+0.15) || 
//						Math.sqrt(x*x+(y+1)*(y+1))<(0.5+0.15) )
//					return 1;
				
				//cylinder4
//				else if(Math.sqrt((x-2.5)*(x-2.5)+(y-0.3)*(y-0.3))<(0.5+0.15) || 
//						Math.sqrt((x+2.5)*(x+2.5)+(y+0.3)*(y+0.3))<(0.5+0.15) )
//					return 1;
				
				//all
				else
					return 0;
			}
		});
		//uv其他边界
		mapNTF_uv.put(NodeType.Neumann, null);
		
		
		HashMap<NodeType, Function> mapNTF_p = new HashMap<NodeType, Function>();
		mapNTF_p.put(NodeType.Dirichlet, new AbstractFunction("x","y") {
			@Override
			public double value(Variable v) {
				double x = v.get("x");
				double y = v.get("y");
				//cylinder1 cylinder2 cylinder4
//				if(Math.abs(x-5.0)<0.01)
//					return 1;
//				else
//					return 0;

				//u_shape
				if(Math.abs(x+5.0)<0.01 && y<0)
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
		
		FEQuadraticV_LinearP fe = new FEQuadraticV_LinearP();
		fe.initDOFIndexGenerator(nodes.size());
		for(int i=1;i<=eList.size();i++) {
			fe.assignTo(eList.at(i));
			//eList.at(i).printDOFInfo();
		}

		//Laplace2D weak form
		WeakFormStokes weakForm = new WeakFormStokes();
		
		//Right hand side(RHS): f = (0,0)'
		weakForm.setF(new SpaceVectorFunction(FC.c0,FC.c0));
		weakForm.setParam(FC.c(1.0),FC.c0);
		//Robin:  k*u_n + d*u - p\mathbf{n} = 0
		VectorFunction d = new SpaceVectorFunction(2);
		d.set(1, FC.c0);
		d.set(2, FC.c0);
		weakForm.setRobin(d);
		
		//Assemble
		AssemblerVector assembler = new AssemblerVector(mesh, weakForm,fe);
		System.out.println("Begin Assemble...");
		assembler.assemble();
		BlockMatrix stiff = (BlockMatrix)assembler.getStiffnessMatrix();
		BlockVector load = (BlockVector)assembler.getLoadVector();
		//stiff.print();
		//load.print();
		//System.out.println(load.norm2());
		//Boundary condition
		VectorFunction diri = new SpaceVectorFunction(3);
		diri.set(1, new AbstractFunction("x","y") {
					@Override
					public double value(Variable v) {
						double x = v.get("x");
						double y = v.get("y");
						
						//cylinder1 cylinder2 cylinder4 
//						if(Math.abs(x+5.0)<0.01)
//							return Math.cos(Math.PI/2*(y/2.5));
						//u_shape
						if(Math.abs(x+5.0)<0.01 && y>0)
							return Math.cos(Math.PI/2*(y-1.5));
						//all
						else
							return 0.0;
					}
				});
		diri.set(2, FC.c0);
		diri.set(3, FC.c0);
		assembler.imposeDirichletCondition(diri);
		load.getBlock(1).print();
		load.getBlock(2).print();
		load.getBlock(3).print();
		System.out.println("Assemble done!");
		
		SchurComplementStokesSolver solver = 
			new SchurComplementStokesSolver(stiff,load);
		
		BlockVector u = solver.solve();
		
		System.out.println("u=");
		for(int i=1;i<=u.getDim();i++)
			System.out.println(String.format("%.3f", u.get(i)));	
	    
		Tools.plotVector(mesh, outputFolder, String.format("%s_uv.dat",file), 
				u.getBlock(1), u.getBlock(2));
		Tools.plotVector(meshOld, outputFolder, String.format("%s_p.dat",file), 
				u.getBlock(3));
		
	}
	
	public static void main(String[] args) {
		cylinder();
	}
}
