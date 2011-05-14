package edu.uta.futureye.test;

import java.util.HashMap;

import edu.uta.futureye.algebra.SolverJBLAS;
import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.core.DOF;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.NodeRefined;
import edu.uta.futureye.core.NodeType;
import edu.uta.futureye.core.Refiner;
import edu.uta.futureye.function.basic.FAxpb;
import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.function.operator.FMath;
import edu.uta.futureye.io.MeshReader;
import edu.uta.futureye.io.MeshWriter;
import edu.uta.futureye.lib.assembler.AssemblerScalar;
import edu.uta.futureye.lib.shapefun.SFBilinearLocal2D;
import edu.uta.futureye.lib.shapefun.SFLinearLocal2D;
import edu.uta.futureye.lib.weakform.WeakFormLaplace2D;
import edu.uta.futureye.util.container.ElementList;

public class TestAdaptive {
	
	
	public static void beforeRefinement() {
//		MeshReader reader = new MeshReader("patch_rectangle.grd");
		MeshReader reader = new MeshReader("patch_rectangle2.grd");
//		MeshReader reader = new MeshReader("patch_rectangle_refine.grd");
		Mesh mesh = reader.read2DMesh();

		mesh.computeNodeBelongsToElements();
		mesh.computeNeighborNodes();
		mesh.computeGlobalEdge();
		mesh.computeNeighborElements();
		
		HashMap<NodeType, Function> mapNTF = new HashMap<NodeType, Function>();
		mapNTF.put(NodeType.Dirichlet, null);	
		mesh.markBorderNode(mapNTF);
		
		SFBilinearLocal2D[] shapeFun = new SFBilinearLocal2D[4];
		for(int i=0;i<4;i++)
			shapeFun[i] = new SFBilinearLocal2D(i+1);
		
		//Asign degree of freedom to element
		for(int i=1;i<=mesh.getElementList().size();i++) {
			Element e = mesh.getElementList().at(i);
			int nDofLocalIndexCounter = 0;
			for(int j=1;j<=e.nodes.size();j++) {
				//Asign shape function to DOF
				DOF dof = new DOF(++nDofLocalIndexCounter,e.nodes.at(j).globalIndex,shapeFun[j-1]);
				e.addNodeDOF(j, dof);
			}
		}
		
		//User defined weak form of PDE (including bounder conditions)
		WeakFormLaplace2D weakForm = new WeakFormLaplace2D();
		//-\Delta{u} = f
		//u(x,y)=0, (x,y)\in\partial{\Omega}
		//\Omega = [0,10]*[0,10]
		//u=[(x-5)^2-25]*[(y-5)^2-25]
		//f=-2*( (x-5)^2 + (y-5)^2 ) + 100
		Function fxm5 = new FAxpb("x",1.0,-5.0);
		Function fym5 = new FAxpb("y",1.0,-5.0);
		weakForm.setF(
				FC.c(-2.0).M(FMath.pow(fxm5, new FC(2.0)) ).A(
						FC.c(-2.0).M(FMath.pow(fym5, new FC(2.0)) )
						).A(FC.c(100.0))
				);
		
		AssemblerScalar assembler = new AssemblerScalar(mesh, weakForm);
		System.out.println("Begin Assemble...");
		assembler.assemble();
		Matrix stiff = assembler.getStiffnessMatrix();
		Vector load = assembler.getLoadVector();
		assembler.imposeDirichletCondition(new FC(0.0));
		System.out.println("Assemble done!");
		
		SolverJBLAS solver = new SolverJBLAS();
		Vector u = solver.solveDGESV(stiff, load);
		
	    System.out.println("u=");
	    for(int i=1;i<=u.getDim();i++)
	        System.out.println(String.format("%.3f", u.get(i)));
	   
	    MeshWriter writer = new MeshWriter(mesh);
	    writer.writeTechplot("patch_rectangle_before_refine.dat", u);			
	}
	
//	public static void patchTest() {
////		MeshReader reader = new MeshReader("patch_rectangle.grd");
//		MeshReader reader = new MeshReader("patch_rectangle2.grd");
////		MeshReader reader = new MeshReader("patch_rectangle_refine.grd");
//		Mesh mesh = reader.read2D();
//		HashMap<NodeType, Function> mapNTF = new HashMap<NodeType, Function>();
//		mapNTF.put(NodeType.Dirichlet, null);	
//
//		mesh.computeNodesBelongToElement();
//		mesh.computeNeiborNode();
//		mesh.markBorderNode(mapNTF);
//		mesh.computeNeighborElement();
//
//		ElementList eList = mesh.getElementList();
//		Element eOld = eList.at(1);
//		ElementList refinedList = eOld.refineOnce();
//		NodeList nList = mesh.getNodeList();
//		
//		//第一步：新增加的结点赋予全局编号，加入mesh对象
//		for(int i=1;i<=refinedList.size();i++) {
//			Element eNew = refinedList.at(i);
//			if(eNew.globalIndex == 0) {
//				eNew.globalIndex = eList.size()+1;
//				eList.add(eNew);
//				for(int j=1;j<=eNew.nodes.size();j++) {
//					Node nNew = eNew.nodes.at(j);
//					if(nNew.globalIndex == 0) {
//						nNew.globalIndex = nList.size()+1;
//						nList.add(nNew);
//					}
//				}
//			}
//		}
//		
//		//第二步：计算hanging node
//		ElementList eNeighbor = eOld.neighbors;
//		for(int i=1;i<=refinedList.size();i++) {
//			Element eRefined = refinedList.at(i);
//			for(int j=1;j<=eRefined.nodes.size();j++) {
//				Node nNew = eRefined.nodes.at(j);
//				if(nNew instanceof NodeRefined) {
//					NodeRefined nRefined = (NodeRefined)nNew;
//					for(int k=1;k<=eNeighbor.size();k++) {
//						EdgeList edges = eNeighbor.at(k).getEdgeList();
//						for(int kk=1;kk<=edges.size();kk++) {
//							if(edges.at(kk).isCoordOnEdge(nRefined.coords())) {
//								NodeList endNodes = edges.at(kk).getEndNodes();
//								nRefined.addConstrainNode(endNodes.at(1));
//								nRefined.addConstrainNode(endNodes.at(2));
//							}
//						}
//					}
//				}
//			}
//		}
//		
//		eList.remove(eOld);
//		
//		mesh.computeNodesBelongToElement();
//		mesh.computeNeiborNode();
//		mesh.computeNeighborElement();
//		mesh.markBorderNode(mapNTF);
//
//		SFBilinearLocal2D[] shapeFun = new SFBilinearLocal2D[4];
//		for(int i=0;i<4;i++)
//			shapeFun[i] = new SFBilinearLocal2D(i+1);
//		SFBilinearLocal2D[] shapeFun2 = new SFBilinearLocal2D[4];
//		for(int i=0;i<4;i++) {
//			shapeFun2[i] = new SFBilinearLocal2D(i+1,0.5);
//		}
//		
//		//Asign degree of freedom to element
//		for(int i=1;i<=mesh.getElementList().size();i++) {
//			Element e = mesh.getElementList().at(i);
//			int nDofLocalIndexCounter = 0;
//			for(int j=1;j<=e.nodes.size();j++) {
//				//Asign shape function to DOF
//				if(e.nodes.at(j) instanceof NodeRefined) {
//					NodeRefined nRefined = (NodeRefined)e.nodes.at(j);
//					if(nRefined.isHangingNode()) {
//						DOF dof  = new DOF(++nDofLocalIndexCounter,nRefined.constrainNodes.at(1).globalIndex,
//								shapeFun2[j-1]);
//						e.addDOF(j, dof);
//						DOF dof2 = new DOF(++nDofLocalIndexCounter,nRefined.constrainNodes.at(2).globalIndex,
//								shapeFun2[j-1]);
//						e.addDOF(j, dof2);
//					} else {
//						DOF dof = new DOF(++nDofLocalIndexCounter,e.nodes.at(j).globalIndex,shapeFun[j-1]);
//						e.addDOF(j, dof);				
//					}
//				} else {
//					DOF dof = new DOF(++nDofLocalIndexCounter,e.nodes.at(j).globalIndex,shapeFun[j-1]);
//					e.addDOF(j, dof);
//				}
//			}
//		}
//		
//		//User defined weak form of PDE (including bounder conditions)
//		WeakFormLaplace2D weakForm = new WeakFormLaplace2D();
//		//-\Delta{u} = f
//		//u(x,y)=0, (x,y)\in\partial{\Omega}
//		//\Omega = [0,10]*[0,10]
//		//u=[(x-5)^2-25]*[(y-5)^2-25]
//		//f=-2*( (x-5)^2 + (y-5)^2 ) + 100
//		Function fxm5 = new FAxpb("x",1.0,-5.0);
//		Function fym5 = new FAxpb("y",1.0,-5.0);
//		weakForm.setF(
//				FOBasic.Plus(
//					FOBasic.Plus(
//						FOBasic.Mult(new FConstant(-2.0), FOBasic.Power(fxm5, new FConstant(2.0)) ),
//						FOBasic.Mult(new FConstant(-2.0), FOBasic.Power(fym5, new FConstant(2.0)) )
//						),new FConstant(100.0)
//					)
//				);
//		
//		Assembler assembler = new Assembler(mesh, weakForm);
//		System.out.println("Begin Assemble...");
//		Matrix stiff = assembler.getStiffnessMatrix();
//		Vector load = assembler.getLoadVector();
//		assembler.imposeDirichletCondition(new FConstant(0.0));
//		System.out.println("Assemble done!");
//		
//		Solver solver = new Solver();
//		Vector u = solver.solve(stiff, load);
//		
//		//hanging node赋值
//		for(int i=1;i<=mesh.getElementList().size();i++) {
//			Element e = mesh.getElementList().at(i);
//			for(int j=1;j<=e.nodes.size();j++) {
//				if(e.nodes.at(j) instanceof NodeRefined) {
//					NodeRefined nRefined = (NodeRefined)e.nodes.at(j);
//					if(nRefined.isHangingNode()) {
//						double hnValue = (u.get(nRefined.constrainNodes.at(1).globalIndex)+
//						u.get(nRefined.constrainNodes.at(2).globalIndex))/2.0;
//						
//						u.set(nRefined.globalIndex, hnValue);
//					}
//				}
//			}
//		}
//		
//	    System.out.println("u=");
//	    for(int i=1;i<=u.getDim();i++)
//	        System.out.println(String.format("%.3f", u.get(i)));
//	   
//	    MeshWriter writer = new MeshWriter(mesh);
//	    writer.writeTechplot("patch_rectangle.dat", u);	
//		
//	}
	
	
	public static void adaptiveTestRectangle() {
//		MeshReader reader = new MeshReader("patch_rectangle.grd");
		MeshReader reader = new MeshReader("patch_rectangle2.grd");
//		MeshReader reader = new MeshReader("patch_rectangle_refine.grd");
	
		Mesh mesh = reader.read2DMesh();
		HashMap<NodeType, Function> mapNTF = new HashMap<NodeType, Function>();
		mapNTF.put(NodeType.Dirichlet, null);	

		mesh.computeNodeBelongsToElements();
		mesh.computeNeighborNodes();
		mesh.markBorderNode(mapNTF);
		mesh.computeGlobalEdge();
		mesh.computeNeighborElements();

		ElementList eList = mesh.getElementList();
		ElementList eToRefine = new ElementList();
		eToRefine.add(eList.at(6));
		eToRefine.add(eList.at(7));
		eToRefine.add(eList.at(10));
		eToRefine.add(eList.at(11));
		
		Refiner.refineOnce(mesh, eToRefine);
		mesh.markBorderNode(mapNTF);
		
		//二次加密
		eToRefine.clear();
		//单元编号的问题该如何处理？
		eToRefine.add(eList.at(17));
		eToRefine.add(eList.at(18));
		//第一步：新增加的结点赋予全局编号，加入mesh对象
		Refiner.refineOnce(mesh, eToRefine);
		mesh.markBorderNode(mapNTF);

		SFBilinearLocal2D[] shapeFun = new SFBilinearLocal2D[4];
		for(int i=0;i<4;i++)
			shapeFun[i] = new SFBilinearLocal2D(i+1);
		SFBilinearLocal2D[] shapeFun2 = new SFBilinearLocal2D[4];
		for(int i=0;i<4;i++) {
			shapeFun2[i] = new SFBilinearLocal2D(i+1,0.5);
		}
		
		//Asign degree of freedom to element
		for(int i=1;i<=mesh.getElementList().size();i++) {
			Element e = mesh.getElementList().at(i);
			int nDofLocalIndexCounter = 0;
			for(int j=1;j<=e.nodes.size();j++) {
				//Asign shape function to DOF
				if(e.nodes.at(j) instanceof NodeRefined) {
					NodeRefined nRefined = (NodeRefined)e.nodes.at(j);
					if(nRefined.isHangingNode()) {
						DOF dof  = new DOF(++nDofLocalIndexCounter,
								nRefined.constrainNodes.at(1).globalIndex,
								shapeFun2[j-1]);
						e.addNodeDOF(j, dof);
						DOF dof2 = new DOF(++nDofLocalIndexCounter,
								nRefined.constrainNodes.at(2).globalIndex,
								shapeFun2[j-1]);
						e.addNodeDOF(j, dof2);
					} else {
						DOF dof = new DOF(++nDofLocalIndexCounter,
								e.nodes.at(j).globalIndex,shapeFun[j-1]);
						e.addNodeDOF(j, dof);				
					}
				} else {
					DOF dof = new DOF(++nDofLocalIndexCounter,
							e.nodes.at(j).globalIndex,shapeFun[j-1]);
					e.addNodeDOF(j, dof);
				}
			}
		}
		
		//User defined weak form of PDE (including bounder conditions)
		WeakFormLaplace2D weakForm = new WeakFormLaplace2D();
		//-\Delta{u} = f
		//u(x,y)=0, (x,y)\in\partial{\Omega}
		//\Omega = [0,10]*[0,10]
		//u=[(x-5)^2-25]*[(y-5)^2-25]
		//f=-2*( (x-5)^2 + (y-5)^2 ) + 100
		Function fxm5 = new FAxpb("x",1.0,-5.0);
		Function fym5 = new FAxpb("y",1.0,-5.0);
		weakForm.setF(
				FC.c(-2.0).M(FMath.pow(fxm5, new FC(2.0)) ).A(
						FC.c(-2.0).M(FMath.pow(fym5, new FC(2.0)) )
						).A(FC.c(100.0))
				);
		
		AssemblerScalar assembler = new AssemblerScalar(mesh, weakForm);
		System.out.println("Begin Assemble...");
		assembler.assemble();
		Matrix stiff = assembler.getStiffnessMatrix();
		Vector load = assembler.getLoadVector();
		assembler.imposeDirichletCondition(new FC(0.0));
		System.out.println("Assemble done!");
		
		SolverJBLAS solver = new SolverJBLAS();
		Vector u = solver.solveDGESV(stiff, load);
		
		//hanging node赋值
		for(int i=1;i<=mesh.getElementList().size();i++) {
			Element e = mesh.getElementList().at(i);
			for(int j=1;j<=e.nodes.size();j++) {
				if(e.nodes.at(j) instanceof NodeRefined) {
					NodeRefined nRefined = (NodeRefined)e.nodes.at(j);
					if(nRefined.isHangingNode()) {
						double hnValue = 
							(u.get(nRefined.constrainNodes.at(1).globalIndex)+
							u.get(nRefined.constrainNodes.at(2).globalIndex))/2.0;
						
						u.set(nRefined.globalIndex, hnValue);
					}
				}
			}
		}
		
	    System.out.println("u=");
	    for(int i=1;i<=u.getDim();i++)
	        System.out.println(String.format("%.3f", u.get(i)));
	   
	    MeshWriter writer = new MeshWriter(mesh);
	    writer.writeTechplot("patch_rectangle.dat", u);	
		
	}

	
	public static void adaptiveTestTriangle() {
		MeshReader reader = new MeshReader("patch_triangle.grd");
	
		Mesh mesh = reader.read2DMesh();
		mesh.computeNodeBelongsToElements();
		mesh.computeNeighborNodes();
		mesh.computeGlobalEdge();
		mesh.computeNeighborElements();
		
		HashMap<NodeType, Function> mapNTF = new HashMap<NodeType, Function>();
		mapNTF.put(NodeType.Dirichlet, null);	
		mesh.markBorderNode(mapNTF);

		ElementList eList = mesh.getElementList();
		ElementList eToRefine = new ElementList();
//直接指定需要加密的单元编号
		eToRefine.add(eList.at(6));
		eToRefine.add(eList.at(16));
		eToRefine.add(eList.at(3));
		eToRefine.add(eList.at(4));
		
		Refiner.refineOnce(mesh, eToRefine);
		mesh.markBorderNode(mapNTF);
		
//		//二次加密
//		eToRefine.clear();
//		//单元编号的问题该如何处理？
//		eToRefine.add(eList.at(17));
//		//第一步：新增加的结点赋予全局编号，加入mesh对象
//		Refiner.refineOnce(mesh, eToRefine);
//		mesh.markBorderNode(mapNTF);

		SFLinearLocal2D[] shapeFun = new SFLinearLocal2D[3];
		for(int i=0;i<3;i++)
			shapeFun[i] = new SFLinearLocal2D(i+1);
		SFLinearLocal2D[] shapeFun2 = new SFLinearLocal2D[3];
		for(int i=0;i<3;i++) {
			shapeFun2[i] = new SFLinearLocal2D(i+1,0.5);
		}
		
		//Asign degree of freedom to element
		for(int i=1;i<=mesh.getElementList().size();i++) {
			Element e = mesh.getElementList().at(i);
			int nDofLocalIndexCounter = 0;
			for(int j=1;j<=e.nodes.size();j++) {
				//Asign shape function to DOF
				if(e.nodes.at(j) instanceof NodeRefined) {
					NodeRefined nRefined = (NodeRefined)e.nodes.at(j);
					if(nRefined.isHangingNode()) {
						DOF dof  = new DOF(++nDofLocalIndexCounter,nRefined.constrainNodes.at(1).globalIndex,
								shapeFun2[j-1]);
						e.addNodeDOF(j, dof);
						DOF dof2 = new DOF(++nDofLocalIndexCounter,nRefined.constrainNodes.at(2).globalIndex,
								shapeFun2[j-1]);
						e.addNodeDOF(j, dof2);
					} else {
						DOF dof = new DOF(++nDofLocalIndexCounter,e.nodes.at(j).globalIndex,shapeFun[j-1]);
						e.addNodeDOF(j, dof);				
					}
				} else {
					DOF dof = new DOF(++nDofLocalIndexCounter,e.nodes.at(j).globalIndex,shapeFun[j-1]);
					e.addNodeDOF(j, dof);
				}
			}
		}
		
		//User defined weak form of PDE (including bounder conditions)
		WeakFormLaplace2D weakForm = new WeakFormLaplace2D();
		//-\Delta{u} = f
		//u(x,y)=0, (x,y)\in\partial{\Omega}
		//\Omega = [0,10]*[0,10]
		//u=[(x-5)^2-25]*[(y-5)^2-25]
		//f=-2*( (x-5)^2 + (y-5)^2 ) + 100
		Function fxm5 = new FAxpb("x",1.0,-5.0);
		Function fym5 = new FAxpb("y",1.0,-5.0);
		weakForm.setF(
						FC.c(-2.0).M(FMath.pow(fxm5, new FC(2.0)) ).A(
						FC.c(-2.0).M(FMath.pow(fym5, new FC(2.0)) )
						).A(FC.c(100.0))
				);
		
		AssemblerScalar assembler = new AssemblerScalar(mesh, weakForm);
		System.out.println("Begin Assemble...");
		assembler.assemble();
		Matrix stiff = assembler.getStiffnessMatrix();
		Vector load = assembler.getLoadVector();
		assembler.imposeDirichletCondition(new FC(0.0));
		System.out.println("Assemble done!");
		
		SolverJBLAS solver = new SolverJBLAS();
		Vector u = solver.solveDGESV(stiff, load);
		
		//hanging node赋值
		for(int i=1;i<=mesh.getElementList().size();i++) {
			Element e = mesh.getElementList().at(i);
			for(int j=1;j<=e.nodes.size();j++) {
				if(e.nodes.at(j) instanceof NodeRefined) {
					NodeRefined nRefined = (NodeRefined)e.nodes.at(j);
					if(nRefined.isHangingNode()) {
						double hnValue = (u.get(nRefined.constrainNodes.at(1).globalIndex)+
						u.get(nRefined.constrainNodes.at(2).globalIndex))/2.0;
						
						u.set(nRefined.globalIndex, hnValue);
					}
				}
			}
		}
		
	    System.out.println("u=");
	    for(int i=1;i<=u.getDim();i++)
	        System.out.println(String.format("%.3f", u.get(i)));
	   
	    MeshWriter writer = new MeshWriter(mesh);
	    writer.writeTechplot("patch_triangle.dat", u);	
		
	}
	
	public static void main(String[] args) {
		beforeRefinement();
		adaptiveTestRectangle();
		adaptiveTestTriangle();
	}
}
