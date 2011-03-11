package edu.uta.futureye.test;

import java.util.HashMap;

import edu.uta.futureye.algebra.Solver;
import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.core.DOF;
import edu.uta.futureye.core.EdgeLocal;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.core.NodeLocal;
import edu.uta.futureye.core.NodeType;
import edu.uta.futureye.core.Vertex;
import edu.uta.futureye.core.intf.Assembler;
import edu.uta.futureye.function.AbstractFunction;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.basic.FX;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.function.operator.FOBasic;
import edu.uta.futureye.io.MeshReader;
import edu.uta.futureye.io.MeshWriter;
import edu.uta.futureye.lib.assembler.AssemblerScalar;
import edu.uta.futureye.lib.assembler.AssemblerScalarFast;
import edu.uta.futureye.lib.element.FEBilinearRectangle;
import edu.uta.futureye.lib.element.FELinearTriangle;
import edu.uta.futureye.lib.shapefun.SFBilinearLocal2D;
import edu.uta.futureye.lib.shapefun.SFLinearLocal2D;
import edu.uta.futureye.lib.shapefun.SFQuadraticLocal2DFast;
import edu.uta.futureye.lib.shapefun.SFSerendipity2D;
import edu.uta.futureye.lib.weakform.WeakFormLaplace2D;
import edu.uta.futureye.util.list.ElementList;
import edu.uta.futureye.util.list.ObjList;

public class LaplaceTest {
	
	public static void triangleTest() {
//		String meshName = "triangle2";
		String meshName = "triangle_refine";
		MeshReader reader = new MeshReader(meshName+".grd");
		Mesh mesh = reader.read2DMesh();
		mesh.computeNodesBelongToElement();
		
		HashMap<NodeType, Function> mapNTF = new HashMap<NodeType, Function>();
//		mapNTF.put(NodeType.Robin, new FAbstract("x","y"){
//			@Override
//			public double value(Variable v) {
//				if(3.0-v.get("x")<0.01)
//					return 1.0;
//				else
//					return -1.0;
//			}
//		});
		mapNTF.put(NodeType.Dirichlet, null);		
		mesh.markBorderNode(mapNTF);

		//Asign degree of freedom to element
/**	
 * Raw Method:
		SFLinearLocal2D[] shapeFun = new SFLinearLocal2D[3];
		for(int i=0;i<3;i++)
			shapeFun[i] = new SFLinearLocal2D(i+1);
		for(int i=1;i<=mesh.getElementList().size();i++) {
			Element e = mesh.getElementList().at(i);
			for(int j=1;j<=e.nodes.size();j++) {
				//Asign shape function to DOF
				DOF dof = new DOF(j,e.nodes.at(j).globalIndex,shapeFun[j-1]);
				e.addNodeDOF(j, dof);
			}
		}
*/		
		//Use element library
		ElementList eList = mesh.getElementList();
		FELinearTriangle linearTriangle = new FELinearTriangle();
		for(int i=1;i<=eList.size();i++)
			linearTriangle.assign(eList.at(i));
		
		//User defined weak form of PDE (including bounder conditions)
		//-\Delta{u} = f
		//u(x,y)=0, (x,y)\in\partial{\Omega}
		//u=(x^2-9)*(y^2-9)
		//f=-2*(x^2+y^2)+36
		WeakFormLaplace2D weakForm = new WeakFormLaplace2D();
		weakForm.setF(
				FOBasic.Plus(
					FOBasic.Plus(
						FOBasic.Mult(new FC(-2.0), FOBasic.Mult(new FX("x"),new FX("x") )),
						FOBasic.Mult(new FC(-2.0), FOBasic.Mult(new FX("y"),new FX("y") ))
						),new FC(36.0)
					)
				);
//		weakForm.setF(new FConstant(-2.0));
		
		weakForm.setParam(
					null,
					null,
					FOBasic.Minus(
							FOBasic.Mult(new FC(6.0), 
							FOBasic.Mult(new FX("y"),new FX("y") )
							),
					new FC(54.0)),null //Robin: 6*y^2-54
				);
		
		Assembler assembler = new AssemblerScalar(mesh, weakForm);
		System.out.println("Begin Assemble...");
		long begin = System.currentTimeMillis();
		assembler.assemble();
		Matrix stiff = assembler.getStiffnessMatrix();
		Vector load = assembler.getLoadVector();
		assembler.imposeDirichletCondition(new FC(0.0));
		long end = System.currentTimeMillis();
		System.out.println("Assemble done!");
		System.out.println("Time used:"+(end-begin));
		
		Solver solver = new Solver();
		Vector u = solver.solve(stiff, load);
	    System.out.println("u=");
	    for(int i=1;i<=u.getDim();i++)
	        System.out.println(String.format("%.3f", u.get(i)));	
	    
	    MeshWriter writer = new MeshWriter(mesh);
	    writer.writeTechplot(meshName+"_out.dat", u);
		
	}
	
	public static void rectangleTest() {
//		MeshReader reader = new MeshReader("rectangle.grd");
		MeshReader reader = new MeshReader("rectangle_refine.grd");
		Mesh mesh = reader.read2DMesh();
		mesh.computeNodesBelongToElement();
		HashMap<NodeType, Function> mapNTF = new HashMap<NodeType, Function>();
//		mapNTF.put(NodeType.Robin, new FAbstract("x","y"){
//			@Override
//			public double value(Variable v) {
//				if(3.0-v.get("x")<0.01)
//					return 1.0;
//				else
//					return -1.0;
//			}
//		});
		mapNTF.put(NodeType.Dirichlet, null);		

		mesh.markBorderNode(mapNTF);
/**
 * Raw Method:
 		SFBilinearLocal2D[] shapeFun = new SFBilinearLocal2D[4];
		shapeFun[0] = new SFBilinearLocal2D(1);
		shapeFun[1] = new SFBilinearLocal2D(2);
		shapeFun[2] = new SFBilinearLocal2D(3);
		shapeFun[3] = new SFBilinearLocal2D(4);
		
		//Asign degree of freedom to element
		for(int i=1;i<=mesh.getElementList().size();i++) {
			Element e = mesh.getElementList().at(i);
			for(int j=1;j<=e.nodes.size();j++) {
				//Asign shape function to DOF
				DOF dof = new DOF(j,e.nodes.at(j).globalIndex,shapeFun[j-1]);
				e.addNodeDOF(j, dof);
			}
		}
*/
		//Use element library
		ElementList eList = mesh.getElementList();
		FEBilinearRectangle bilinearRectangle = new FEBilinearRectangle();
		for(int i=1;i<=eList.size();i++)
			bilinearRectangle.assign(eList.at(i));
		
		//User defined weak form of PDE (including bounder conditions)
		WeakFormLaplace2D weakForm = new WeakFormLaplace2D();
		weakForm.setF(
				FOBasic.Plus(
					FOBasic.Plus(
						FOBasic.Mult(new FC(-2.0), FOBasic.Mult(new FX("x"),new FX("x") )),
						FOBasic.Mult(new FC(-2.0), FOBasic.Mult(new FX("y"),new FX("y") ))
						),new FC(36.0)
					)
				);
		
		weakForm.setParam(
				null,
				null,
				new FC(0.05),null //Robin: d*u + k*u_n = q
				); 	
		
		Assembler assembler = new AssemblerScalar(mesh, weakForm);
		System.out.println("Begin Assemble...");
		assembler.assemble();
		Matrix stiff = assembler.getStiffnessMatrix();
		Vector load = assembler.getLoadVector();
		assembler.imposeDirichletCondition(new FC(0.0));
		System.out.println("Assemble done!");
		
		Solver solver = new Solver();
		Vector u = solver.solve(stiff, load);
	    System.out.println("u=");
	    for(int i=1;i<=u.getDim();i++)
	        System.out.println(String.format("%.3f", u.get(i)));
	   
	    MeshWriter writer = new MeshWriter(mesh);
	    writer.writeTechplot("rectangle_refine2_out.dat", u);
	}
	
	
	public static void mixedTest() {
		//注：
		//数值积分三角形上选择2阶，四边形上选择1阶
		//三角形上选择4阶，四边形上选择2阶（bugfix）
		MeshReader reader = new MeshReader("mixed.grd");
		Mesh mesh = reader.read2DMesh();
		mesh.computeNodesBelongToElement();
		HashMap<NodeType, Function> mapNTF = new HashMap<NodeType, Function>();
//		mapNTF.put(NodeType.Robin, new FAbstract("x","y"){
//			@Override
//			public double value(Variable v) {
//				if(3.0-v.get("x")<0.01)
//					return 1.0;
//				else
//					return -1.0;
//			}
//		});
		mapNTF.put(NodeType.Dirichlet, null);		

		mesh.markBorderNode(mapNTF);

		SFBilinearLocal2D[] shapeFunRect = new SFBilinearLocal2D[4];
		shapeFunRect[0] = new SFBilinearLocal2D(1);
		shapeFunRect[1] = new SFBilinearLocal2D(2);
		shapeFunRect[2] = new SFBilinearLocal2D(3);
		shapeFunRect[3] = new SFBilinearLocal2D(4);
		
		SFLinearLocal2D[] shapeFunTri = new SFLinearLocal2D[3];
		shapeFunTri[0] = new SFLinearLocal2D(1);
		shapeFunTri[1] = new SFLinearLocal2D(2);
		shapeFunTri[2] = new SFLinearLocal2D(3);
		
		//Asign degree of freedom to element
		for(int i=1;i<=mesh.getElementList().size();i++) {
			Element e = mesh.getElementList().at(i);
			if(e.nodes.size() == 4) {
				for(int j=1;j<=e.nodes.size();j++) {
					//Asign shape function to DOF
					DOF dof = new DOF(j,e.nodes.at(j).globalIndex,shapeFunRect[j-1]);
					e.addNodeDOF(j, dof);
				}
			} else if(e.nodes.size() == 3) {
				for(int j=1;j<=e.nodes.size();j++) {
					//Asign shape function to DOF
					DOF dof = new DOF(j,e.nodes.at(j).globalIndex,shapeFunTri[j-1]);
					e.addNodeDOF(j, dof);
				}
			} else {
				System.out.println("Error: e.nodes.size()="+e.nodes.size());
			}
		}
		
		//User defined weak form of PDE (including bounder conditions)
		WeakFormLaplace2D weakForm = new WeakFormLaplace2D();
		weakForm.setF(
				FOBasic.Plus(
					FOBasic.Plus(
						FOBasic.Mult(new FC(-2.0), FOBasic.Mult(new FX("x"),new FX("x") )),
						FOBasic.Mult(new FC(-2.0), FOBasic.Mult(new FX("y"),new FX("y") ))
						),new FC(36.0)
					)
				);
		weakForm.setParam(
				null,
				null,
				new FC(0.05),null //Robin: d*u + k*u_n = q
				);
		Assembler assembler = new AssemblerScalarFast(mesh, weakForm);
		System.out.println("Begin Assemble...");
		assembler.assemble();
		Matrix stiff = assembler.getStiffnessMatrix();
		Vector load = assembler.getLoadVector();
		assembler.imposeDirichletCondition(new FC(0.0));
		System.out.println("Assemble done!");
		
		Solver solver = new Solver();
		Vector u = solver.solve(stiff, load);
	    System.out.println("u=");
	    for(int i=1;i<=u.getDim();i++)
	        System.out.println(String.format("%.3f", u.get(i)));
	    
	    MeshWriter writer = new MeshWriter(mesh);
	    writer.writeTechplot("mixed_out.dat", u);
	    
	}
	
	
	public static void serendipityTest() {
		MeshReader reader = new MeshReader("rectangle.grd");
//		MeshReader reader = new MeshReader("rectangle_refine.grd");
		Mesh mesh = reader.read2DMesh();
		
		//Add nodes for serendipity element
		//int indexSet[] = {1,2,3,4,1};

		for(int i=1;i<=mesh.getElementList().size();i++) {
			Element e = mesh.getElementList().at(i);
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
		
//			for(int j=1;j<=4;j++) {
//				Node l = e.nodes.at(indexSet[j-1]);
//				Node r = e.nodes.at(indexSet[j]);
//				double cx = (l.coord(1)+r.coord(1))/2.0;
//				double cy = (l.coord(2)+r.coord(2))/2.0;
//				Node node = new Node(2);
//				node.set(mesh.getNodeList().size()+1, cx,cy );
//				Node findNode = mesh.containNode(node);
//				if(findNode == null) {
//					e.addVertexNode(node, false);
//					mesh.addNode(node);
//				} else {
//					e.addVertexNode(findNode, false);
//				}
//			}
		}
		
		mesh.computeNodesBelongToElement();
		HashMap<NodeType, Function> mapNTF = new HashMap<NodeType, Function>();
		mapNTF.put(NodeType.Dirichlet, null);
		mesh.markBorderNode(mapNTF);

		
		SFSerendipity2D[] shapeFun = new SFSerendipity2D[8];
		for(int i=1;i<=8;i++)
			shapeFun[i-1] = new SFSerendipity2D(i);
		
		//Asign degree of freedom to element
		for(int i=1;i<=mesh.getElementList().size();i++) {
			Element e = mesh.getElementList().at(i);
			for(int j=1;j<=e.nodes.size();j++) {
				//Asign shape function to DOF
				DOF dof = new DOF(j,e.nodes.at(j).globalIndex,shapeFun[j-1]);
				e.addNodeDOF(j, dof);
			}
		}
		
		//User defined weak form of PDE (including bounder conditions)
		WeakFormLaplace2D weakForm = new WeakFormLaplace2D();
		weakForm.setF(
				FOBasic.Plus(
					FOBasic.Plus(
						FOBasic.Mult(new FC(-2.0), FOBasic.Mult(new FX("x"),new FX("x") )),
						FOBasic.Mult(new FC(-2.0), FOBasic.Mult(new FX("y"),new FX("y") ))
						),new FC(36.0)
					)
				);
		
		weakForm.setParam(
				null,
				null,
				new FC(0.05),null //Robin: d*u + k*u_n = q
				); 
		
		Assembler assembler = new AssemblerScalar(mesh, weakForm);
		System.out.println("Begin Assemble...");
		assembler.assemble();
		Matrix stiff = assembler.getStiffnessMatrix();
		Vector load = assembler.getLoadVector();
		assembler.imposeDirichletCondition(new FC(0.0));
		System.out.println("Assemble done!");
		
		Solver solver = new Solver();
		Vector u = solver.solve(stiff, load);
	    System.out.println("u=");
	    for(int i=1;i<=u.getDim();i++)
	        System.out.println(String.format("%.3f", u.get(i)));
	    
	    MeshWriter writer = new MeshWriter(mesh);
	    writer.writeTechplot("serendipity.dat", u);	    
	}
	
	public static void quadraticLocal2DTest() {
		MeshReader reader = new MeshReader("triangle.grd");
//		MeshReader reader = new MeshReader("triangle_refine0.grd");
//		MeshReader reader = new MeshReader("triangle_refine2.grd");
		Mesh mesh = reader.read2DMesh();
		
		//Add nodes for quadratic element
		//int indexSet[] = {1,2,3,1};
		for(int i=1;i<=mesh.getElementList().size();i++) {
			Element e = mesh.getElementList().at(i);
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
			
			
//			for(int j=1;j<=3;j++) {
//				Node l = e.nodes.at(indexSet[j-1]);
//				Node r = e.nodes.at(indexSet[j]);
//				double cx = (l.coord(1)+r.coord(1))/2.0;
//				double cy = (l.coord(2)+r.coord(2))/2.0;
//				Node node = new Node(2);
//				node.set(mesh.getNodeList().size()+1, cx,cy );
//				Node findNode = mesh.containNode(node);
//				if(findNode == null) {
//					e.addVertexNode(node, false);
//					mesh.addNode(node);
//				} else {
//					e.addVertexNode(findNode, false);
//				}
//			}
		}
		
		mesh.computeNodesBelongToElement();
		HashMap<NodeType, Function> mapNTF = new HashMap<NodeType, Function>();
//		mapNTF.put(NodeType.Robin, new FAbstract("x","y"){
//									@Override
//									public double value(Variable v) {
//										if(3.0-v.get("x")<0.01)
//											return 1.0;
//										else
//											return -1.0;
//									}
//									});
		mapNTF.put(NodeType.Dirichlet, null);
		mesh.markBorderNode(mapNTF);

	
//		SFLinearLocal2D[] shapeFun = new SFLinearLocal2D[3];
//		for(int i=0;i<3;i++)
//			shapeFun[i] = new SFLinearLocal2D(i+1);
		SFQuadraticLocal2DFast[] shapeFun = new SFQuadraticLocal2DFast[6];
		for(int i=1;i<=6;i++)
			shapeFun[i-1] = new SFQuadraticLocal2DFast(i);
//		SFQuadraticLocal2D[] shapeFun = new SFQuadraticLocal2D[6];
//		for(int i=1;i<=6;i++)
//			shapeFun[i-1] = new SFQuadraticLocal2D(i);
		
		//Asign degree of freedom to element
		for(int i=1;i<=mesh.getElementList().size();i++) {
			Element e = mesh.getElementList().at(i);
			for(int j=1;j<=e.nodes.size();j++) {
				//Asign shape function to DOF
				DOF dof = new DOF(j,e.nodes.at(j).globalIndex,shapeFun[j-1]);
				e.addNodeDOF(j, dof);
			}
		}
		
		//User defined weak form of PDE (including bounder conditions)
		WeakFormLaplace2D weakForm = new WeakFormLaplace2D();
		weakForm.setF(
				FOBasic.Plus(
					FOBasic.Plus(
						FOBasic.Mult(new FC(-2.0), FOBasic.Mult(new FX("x"),new FX("x") )),
						FOBasic.Mult(new FC(-2.0), FOBasic.Mult(new FX("y"),new FX("y") ))
						),new FC(36.0)
					)
				);
		
		weakForm.setParam(
				null,
				null,
				FOBasic.Minus(
						FOBasic.Mult(new FC(6.0), 
						FOBasic.Mult(new FX("y"),new FX("y") )
						),
				new FC(54.0)),null //Robin: 6*y^2-54
				); 
		
		Assembler assembler = new AssemblerScalar(mesh, weakForm);
		System.out.println("Begin Assemble...");
		assembler.assemble();
		Matrix stiff = assembler.getStiffnessMatrix();
		Vector load = assembler.getLoadVector();
		assembler.imposeDirichletCondition(new FC(0.0));
		System.out.println("Assemble done!");
		
		Solver solver = new Solver();
		Vector u = solver.solve(stiff, load);
	    System.out.println("u=");
	    for(int i=1;i<=u.getDim();i++)
	        System.out.println(String.format("%.3f", u.get(i)));
	    
	    MeshWriter writer = new MeshWriter(mesh);
	    writer.writeTechplot("quadraticLocal2D.dat", u);	    
	}
	
	
	/**
	 * Border Depth Compensation test
	 * 
	 * Dirichlet边界条件其中一个结点上的值比其他结点都小，
	 * 结算其影响范围
	 */
	public static void triangleBDCTest() {
		//String meshName = "triangle_refine";
		String meshName = "triangle";
		MeshReader reader = new MeshReader(meshName+".grd");
		Mesh mesh = reader.read2DMesh();
		mesh.computeNodesBelongToElement();
		
		HashMap<NodeType, Function> mapNTF = new HashMap<NodeType, Function>();
		mapNTF.put(NodeType.Dirichlet, null);		
		mesh.markBorderNode(mapNTF);

		SFLinearLocal2D[] shapeFun = new SFLinearLocal2D[3];
		for(int i=0;i<3;i++)
			shapeFun[i] = new SFLinearLocal2D(i+1);
		
		//Asign degree of freedom to element
		for(int i=1;i<=mesh.getElementList().size();i++) {
			Element e = mesh.getElementList().at(i);
			for(int j=1;j<=e.nodes.size();j++) {
				//Asign shape function to DOF
				DOF dof = new DOF(j,e.nodes.at(j).globalIndex,shapeFun[j-1]);
				e.addNodeDOF(j, dof);
			}
		}
		
		//User defined weak form of PDE (including bounder conditions)
		//-\Delta{u} = f
		//u(x,y)=0, (x,y)\in\partial{\Omega}
		//u=(x^2-9)*(y^2-9)
		//f=-2*(x^2+y^2)+36
		WeakFormLaplace2D weakForm = new WeakFormLaplace2D();
		weakForm.setF(
				FOBasic.Plus(
					FOBasic.Plus(
						FOBasic.Mult(new FC(-2.0), FOBasic.Mult(new FX("x"),new FX("x") )),
						FOBasic.Mult(new FC(-2.0), FOBasic.Mult(new FX("y"),new FX("y") ))
						),new FC(36.0)
					)
				);
		
		Assembler assembler = new AssemblerScalar(mesh, weakForm);
		System.out.println("Begin Assemble...");
		long begin = System.currentTimeMillis();
		assembler.assemble();
		Matrix stiff = assembler.getStiffnessMatrix();
		Vector load = assembler.getLoadVector();
		assembler.imposeDirichletCondition(new AbstractFunction("x","y"){
			@Override
			public double value(Variable v) {
				if(3.0-v.get("x")<0.01 && 
						Math.abs(1.0-v.get("y"))<0.2)
					return 0.0;
				else
					return 10.0;
			}
		});
		long end = System.currentTimeMillis();
		System.out.println("Assemble done!");
		System.out.println("Time used:"+(end-begin));
		
		Solver solver = new Solver();
		Vector u = solver.solve(stiff, load);
	    System.out.println("u=");
	    for(int i=1;i<=u.getDim();i++)
	        System.out.println(String.format("%.3f", u.get(i)));	
	    
	    MeshWriter writer = new MeshWriter(mesh);
	    writer.writeTechplot(meshName+"_BDC_out.dat", u);
		
	}
	public static void main(String[] args) {
//		triangleTest();
//		rectangleTest();
	
//TODO
//		mixedTest(); //OK 数值积分的问题fixed
		
//		serendipityTest();
		quadraticLocal2DTest();
//		triangleBDCTest();
	}
		
}