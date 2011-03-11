package edu.uta.futureye.prostate;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;

import edu.uta.futureye.algebra.Solver;
import edu.uta.futureye.algebra.SparseVector;
import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.core.DOF;
import edu.uta.futureye.core.EdgeLocal;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.core.NodeLocal;
import edu.uta.futureye.core.NodeRefined;
import edu.uta.futureye.core.NodeType;
import edu.uta.futureye.core.Refiner;
import edu.uta.futureye.core.Vertex;
import edu.uta.futureye.core.intf.Assembler;
import edu.uta.futureye.function.AbstractFunction;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.basic.FDelta;
import edu.uta.futureye.function.basic.Vector2Function;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.function.intf.ScalarShapeFunction;
import edu.uta.futureye.function.operator.FMath;
import edu.uta.futureye.function.operator.FOBasic;
import edu.uta.futureye.io.MeshReader;
import edu.uta.futureye.io.MeshWriter;
import edu.uta.futureye.lib.assembler.AssemblerScalar;
import edu.uta.futureye.lib.assembler.AssemblerScalarFast;
import edu.uta.futureye.lib.shapefun.SFBilinearLocal2D;
import edu.uta.futureye.lib.shapefun.SFLinearLocal2D;
import edu.uta.futureye.lib.shapefun.SFQuadraticLocal2D;
import edu.uta.futureye.lib.weakform.WeakFormDerivative;
import edu.uta.futureye.lib.weakform.WeakFormL22D;
import edu.uta.futureye.lib.weakform.WeakFormLaplace2D;
import edu.uta.futureye.util.Constant;
import edu.uta.futureye.util.PairDoubleInteger;
import edu.uta.futureye.util.Utils;
import edu.uta.futureye.util.list.ElementList;
import edu.uta.futureye.util.list.NodeList;
import edu.uta.futureye.util.list.ObjList;

/**
 * Solve: 
 * 1.正问题(Forward problem):
 *   Laplace(u) - a*u = -\delta, in \Omega
 *   u = u0,                     on \Gamma1
 *   u + u_n = 0,                on \Gamma2
 *=> Multiply by '-k'
 *   -k*Laplace(u) + c*u = k*\delta
 *=>
 *   A(u, v) = (f, v)
 * where
 *   A(u, v) = ((k*u_x, v_x) + (k*u_y, v_y) ) - (k*u_n,v)_\Gamma2 + (c*u, v)
 *=>
 *   A(u, v) = ((k*u_x, v_x) + (k*u_y, v_y) ) + (k*u,v)_\Gamma2 + (c*u, v)
 * where
 *   \Gamma1: Dirichlet boundary of \Omega
 *   \Gamma2: Neumann(Robin) boundary of \Omega
 *   u_n: \frac{\pratial{u}}{\partial{n}}
 *   n: unit norm vector of \Omega
 * Parameters:
 *   a(x,y) = 3*mu_a*mu_s'
 *   k(x,y) = 1/(3*mu_s')
 *   c(x,y) = mu_a
 *   \delta(x,y) = delta function, light source
 *
 * 2.系数反问题(Parameter inverse problem)
 *   (U*u, v) = (f, v) - (k*grad(U),grad(v))
 * where 
 *   u=mu_a is unknown
 * Parameters:
 *   k(x,y): 1/(3*mu_s')
 *   U(x,y): solution of forward problem
 *   f(x,y): delta function, light source
 *   
 * @author liuyueming
 */
public class Model {
	public String outputFolder = "";
	
	//Light source
	public Function delta = null;
	public Variable lightSource = null; //light source position
	
	//Inclusion mu_a
	public Function mu_a = null;
	
	//Inclusion 1/(3*mu_s') = 1.0/30.0 ?
	public Function k = new FC(0.02);
	
	/**
	 * type=1: one inclusion
	 * type=2: two inclusion
	 * @param incX
	 * @param incY
	 * @param incR
	 * @param maxMu_a
	 * @param type
	 */
	public void setMu_a(double incX, double incY, double incR, double maxMu_a,
			int type) {
		final double fcx = incX;
		final double fcy = incY;
		final double fcr = incR;
		final double fmu_a = maxMu_a;
		final double distance = 0.8;
		if(type == 1) {
			mu_a = new AbstractFunction("x","y"){
				@Override
				public double value(Variable v) {
					double bk = 0.1;
					double dx = v.get("x")-fcx;
					double dy = v.get("y")-fcy;
					if(Math.sqrt(dx*dx+dy*dy) < fcr) {
						double r = fmu_a*Math.cos((Math.PI/2)*Math.sqrt(dx*dx+dy*dy)/fcr); 
						return r<bk?bk:r;
					}
					else
						return bk;
				}
			};
		} else if(type == 2) {
			mu_a = new AbstractFunction("x","y"){
				@Override
				public double value(Variable v) {
					double bk = 0.1;
					double dx = v.get("x")-fcx;
					double dy = v.get("y")-fcy;
					double dx1 = v.get("x")-(fcx+distance);
					double dy1 = v.get("y")-fcy;
					if(Math.sqrt(dx*dx+dy*dy) < fcr) {
						double r = fmu_a*Math.cos((Math.PI/2)*Math.sqrt(dx*dx+dy*dy)/fcr); 
						return r<bk?bk:r;
					}
					else if(Math.sqrt(dx1*dx1+dy1*dy1) < fcr) {
						double r = fmu_a*Math.cos((Math.PI/2)*Math.sqrt(dx1*dx1+dy1*dy1)/fcr); 
						return r<bk?bk:r;
					}
					else
						return bk;
				}
			};			
		}
	}
	
	public void setDelta(double x,double y) {
		this.lightSource = new Variable();
		this.lightSource.set("x", x);
		this.lightSource.set("y", y);
		//delta = new FDelta(this.lightSource,0.01,2e5);
		//测试将dleta函数变得平缓
		delta = new FDelta(this.lightSource,0.05,2e5);
	}

	public Vector solveForwardNeumann(Mesh mesh) {
		//Mark border type
		HashMap<NodeType, Function> mapNTF = new HashMap<NodeType, Function>();
		mapNTF.clear();
		mapNTF.put(NodeType.Robin, null);
		mesh.clearBorderNodeMark();
		mesh.markBorderNode(mapNTF);
		
		WeakFormLaplace2D weakForm = new WeakFormLaplace2D();
		
		//Right hand side
		weakForm.setF(FMath.Mult(this.k,this.delta));
		
		// *** u + u_n = 0, on \Gamma2 ***
		//   A(u, v) = ((k*u_x, v_x) + (k*u_y, v_y) ) - (k*u_n,v)_\Gamma2 + (c*u, v)
		//=>
		//   A(u, v) = ((k*u_x, v_x) + (k*u_y, v_y) ) + (k*u,v)_\Gamma2 + (c*u, v)
		//Robin: d*u + k*u_n = q
		weakForm.setParam(
				//ERROR: this.k, this.mu_a, new FConstant(0.01),null //first test
				//ERROR: this.k, this.mu_a, null, null //k*u_n=0
				//ERROR: this.k, this.mu_a, null, this.mu_a //d==c, q=0
				//ERROR: this.k, this.mu_a, FMath.Mult(this.k,this.delta), null //k*u_n = q, q=f
				this.k, this.mu_a, null, this.k //d==k,q=0 (即：u_n + u =0)
			);
		
		Assembler assembler = new AssemblerScalar(mesh, weakForm);
		System.out.println("Begin Assemble...solveForwardNeumann");
		assembler.assemble();
		Matrix stiff = assembler.getStiffnessMatrix();
		Vector load = assembler.getLoadVector();
		assembler.imposeDirichletCondition(new FC(0.0));
		System.out.println("Assemble done!");

		Solver solver = new Solver();
		Vector u = solver.solve(stiff, load);
		return u;
	}	

	public Vector solveForwardDirichlet(Mesh mesh, Function diri) {
		//Mark border type
		HashMap<NodeType, Function> mapNTF = new HashMap<NodeType, Function>();
		mapNTF.clear();
		
		mapNTF.put(NodeType.Dirichlet, new AbstractFunction("x","y"){
			@Override
			public double value(Variable v) {
				//double x = v.get("x");
				double y = v.get("y");
				if(Math.abs(y - 3.0) < Constant.eps)
					return 1.0;
				else
					return -1.0;
			}
		});
		
		mapNTF.put(NodeType.Robin, null);
		mesh.clearBorderNodeMark();
		mesh.markBorderNode(mapNTF);
		
		WeakFormLaplace2D weakForm = new WeakFormLaplace2D();
		
		//Right hand side
		weakForm.setF(FMath.Mult(this.k,this.delta));

		weakForm.setParam(
				this.k, this.mu_a, null, this.k //d==k,q=0 (即：u_n + u =0)
			);
		
		Assembler assembler = new AssemblerScalarFast(mesh, weakForm);
		System.out.println("Begin Assemble...solveForwardDirichlet");
		assembler.assemble();
		Matrix stiff = assembler.getStiffnessMatrix();
		Vector load = assembler.getLoadVector();
		//Dirichlet condition
		assembler.imposeDirichletCondition(diri);
		
		System.out.println("Assemble done!");

		Solver solver = new Solver();
		Vector u = solver.solve(stiff, load);
		return u;
	}	
	
	public Vector solveParamInverse(Mesh mesh, Vector U) {
		HashMap<NodeType, Function> mapNTF2 = new HashMap<NodeType, Function>();
		mapNTF2.put(NodeType.Dirichlet, null);
		mesh.clearBorderNodeMark();
		mesh.markBorderNode(mapNTF2);
		
		//Weak form
		WeakFormL22D weakFormL2 = new WeakFormL22D();
		//Right hand side
		weakFormL2.setF(this.delta);
		//Parameters
		weakFormL2.setParam(
				this.k, new Vector2Function(U)
			);
		
		Assembler assembler = new AssemblerScalar(mesh, weakFormL2);
		System.out.println("Begin Assemble...solveParamInverse");
		assembler.assemble();
		Matrix stiff = assembler.getStiffnessMatrix();
		Vector load = assembler.getLoadVector();
		assembler.imposeDirichletCondition(new FC(0.1));
		System.out.println("Assemble done!");
		
		Solver solver = new Solver();
		Vector u = solver.solve(stiff, load);
		return u;
	}
	
	
	public Vector solveEnhance(Mesh mesh, 
			Vector a_m, Vector a_m_1, Vector u_m_1,
			int iterNum
			) {
		//Mark border type
		HashMap<NodeType, Function> mapNTF = new HashMap<NodeType, Function>();
		mapNTF.clear();
		mapNTF.put(NodeType.Dirichlet, null);
		mesh.clearBorderNodeMark();
		mesh.markBorderNode(mapNTF);
		
		WeakFormLaplace2D weakForm = new WeakFormLaplace2D();
		
		//Right hand side
		Function lamd = 
				FOBasic.Mult(
						new FC(Math.PI*Math.PI*Math.pow(Math.E, iterNum-1)),
						FOBasic.Power(
								FOBasic.Minus(new Vector2Function(a_m), new Vector2Function(a_m_1)),
								new FC(2.0)
						)
				);
		//lamd = new FConstant(0.8);
		
		Function rhs = 
			FOBasic.Mult(
				FOBasic.Mult(FOBasic.Power(new FC(Math.E), lamd),new FC(Math.pow(1.05, -iterNum))),
				FOBasic.Mult(
					FOBasic.Minus(new Vector2Function(a_m), new Vector2Function(a_m_1)),
					new Vector2Function(u_m_1)
				));
		
		this.plotFunction(mesh, rhs, "prostate_test1_enhance_rhs"+iterNum+".dat");
		weakForm.setF(rhs);
		
		weakForm.setParam(
			new FC(1.0),
			FOBasic.Minus(new FC(0.0), new Vector2Function(a_m)),
			null,null
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
		return u;
	}

	public void plotFunction(Mesh mesh, Function fun, String fileName) {
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
	
	public Vector discreteFunction(Mesh mesh,Function fun) {
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
	    return v;
	}
	
	public void plotVector(Mesh mesh, Vector v, String fileName) {
	    MeshWriter writer = new MeshWriter(mesh);
	    if(!outputFolder.isEmpty()) {
		    File file = new File("./"+outputFolder);
			if(!file.exists()) {
				file.mkdirs();
			}
	    }
	    writer.writeTechplot("./"+outputFolder+"/"+fileName, v);
	}
	
	/**
	 * 依赖于mesh的大小！！！
	 * 依赖于网格边界结点的边界条件标记！！！
	 * @param mesh
	 * @param bkUL
	 * @param incUL
	 * @return
	 */
	public Vector computeTailLeftLightSource(Mesh mesh,Vector bkUL, Vector incUL) {
	    double lightX = this.lightSource.get("x");
	    double lightY = this.lightSource.get("y");
	    NodeList list = mesh.getNodeList();
	    int nNode = list.size();
	    Vector rlt = new SparseVector(nNode);
	    for(int i=1;i<=nNode;i++) {
	    	Node node = list.at(i);
	    	rlt.set(i, bkUL.get(i));
	    	if(node.coord(1) > lightX && node.coord(1) < 5.0-lightX) {
	    		double dis = Math.sqrt( (lightX-node.coord(1))*(lightX-node.coord(1)) +
	    				(lightY-node.coord(2))*(lightY-node.coord(2)) );
	    		double[] coord = {lightX+dis, 3.0};
	    		if(lightX+dis < 5.0) {
		    		Element e = mesh.getElementByCoord(coord);
		    		if(e != null) {
			    		NodeList borderNodes = e.getNodesByType(NodeType.Robin);
			    		if(borderNodes.size()>0) {
			    			double r = 0.0;
			    			//if(borderNodes.size() == 2 ) {
				    			Node p1 = borderNodes.at(1);
				    			Node p2 = borderNodes.at(2);
				    			Node p = new Node(p1.dim());
				    			p.set(0, coord);
				    			r = Utils.linearInterpolate(
				    					p1,p2,p,
				    					incUL.get(p1.globalIndex), incUL.get(p2.globalIndex));
				    			
//				    		//TODO 用结点个数判断是否二次插值，不靠谱！！！
//			    			} else if(borderNodes.size() == 3){
//			    				Node p1 = borderNodes.at(1);
//				    			Node p2 = borderNodes.at(3); //!!!
//				    			Node p3 = borderNodes.at(2);
//				    			Node p = new Node(p1.dim());
//				    			p.set(0, coord);
//				    			r = Utils.quadraticInterpolate(
//				    					p1,p2,p3,p,
//				    					incUL.get(p1.globalIndex), incUL.get(p2.globalIndex),
//				    					incUL.get(p3.globalIndex));
//			    			}
			    			rlt.set(i, r);
			    		}
		    		}
	    		}
	    	}
	    }
	    return rlt;	
	}

	
	public Vector computeTailLeftLightSource2(Mesh mesh,Vector bkUL, Vector incUL) {
	    double lightX = this.lightSource.get("x");
	    double lightY = this.lightSource.get("y");
	    NodeList list = mesh.getNodeList();
	    int nNode = list.size();
	    Vector rlt = new SparseVector(nNode);
	    for(int i=1;i<=nNode;i++) {
	    	Node node = list.at(i);
	    	rlt.set(i, bkUL.get(i));
	    	if(node.coord(1) > lightX && node.coord(1) < 5.0-lightX) {
	    		double dis = Math.sqrt( (lightX-node.coord(1))*(lightX-node.coord(1)) +
	    				(lightY-node.coord(2))*(lightY-node.coord(2)) );
	    		double[] coord = {lightX+dis, 3.0};
	    		if(lightX+dis < 5.0) {
		    		Element e = mesh.getElementByCoord(coord);
		    		if(e != null) {
		    			//用顶点的坐标来判断，不再使用边界条件。因为自适应加密网格后，如果边界单元被加密
		    			//边界上的边上的加密结点是按照内点条件处理的。
			    		NodeList borderNodes = e.getNodesByType(NodeType.Robin);
			    		if(borderNodes.size()>0) {
			    			double r = 0.0;
			    			//if(borderNodes.size() == 2 ) {
				    			Node p1 = borderNodes.at(1);
				    			Node p2 = borderNodes.at(2);
				    			Node p = new Node(p1.dim());
				    			p.set(0, coord);
				    			r = Utils.linearInterpolate(
				    					p1,p2,p,
				    					incUL.get(p1.globalIndex), incUL.get(p2.globalIndex));
				    			
//				    		//TODO 用结点个数判断是否二次插值，不靠谱！！！
//			    			} else if(borderNodes.size() == 3){
//			    				Node p1 = borderNodes.at(1);
//				    			Node p2 = borderNodes.at(3); //!!!
//				    			Node p3 = borderNodes.at(2);
//				    			Node p = new Node(p1.dim());
//				    			p.set(0, coord);
//				    			r = Utils.quadraticInterpolate(
//				    					p1,p2,p3,p,
//				    					incUL.get(p1.globalIndex), incUL.get(p2.globalIndex),
//				    					incUL.get(p3.globalIndex));
//			    			}
			    			rlt.set(i, r);
			    		}
		    		}
	    		}
	    	}
	    }
	    return rlt;	
	}
	
	/**
	 * 依赖于mesh的大小！！！
	 * 依赖于网格边界结点的边界条件标记！！！
	 * @param mesh
	 * @param bkUR
	 * @param incUR
	 * @return
	 */
	public Vector computeTailRightLightSource(Mesh mesh,Vector bkUR, Vector incUR) {
	    double lightX = this.lightSource.get("x");
	    double lightY = this.lightSource.get("y");
	    NodeList list = mesh.getNodeList();
	    int nNode = list.size();
	    Vector rlt = new SparseVector(nNode);
	    for(int i=1;i<=nNode;i++) {
	    	Node node = list.at(i);
	    	rlt.set(i, bkUR.get(i));
	    	if(node.coord(1) < lightX && node.coord(1) > 5.0-lightX) {
	    		double dis = Math.sqrt( (lightX-node.coord(1))*(lightX-node.coord(1)) +
	    				(lightY-node.coord(2))*(lightY-node.coord(2)) );
	    		double[] coord = {lightX-dis, 3.0};
	    		if(lightX-dis > 0.0) {
		    		Element e = mesh.getElementByCoord(coord);
		    		if(e != null) {
			    		NodeList borderNodes = e.getNodesByType(NodeType.Robin);
			    		if(borderNodes.size()>0) {
			    			double r = 0.0;
			    			//if(borderNodes.size() == 2) {
				    			Node p1 = borderNodes.at(1);
				    			Node p2 = borderNodes.at(2);
				    			Node p = new Node(p1.dim());
				    			p.set(0, coord);
				    			r = Utils.linearInterpolate(
				    					p1,p2,p,
				    					incUR.get(p1.globalIndex), incUR.get(p2.globalIndex));
//			    			} else if(borderNodes.size() == 3){
//			    				Node p1 = borderNodes.at(1);
//				    			Node p2 = borderNodes.at(3); //!!!
//				    			Node p3 = borderNodes.at(2);
//				    			Node p = new Node(p1.dim());
//				    			p.set(0, coord);
//				    			r = Utils.quadraticInterpolate(
//				    					p1,p2,p3,p,
//				    					incUR.get(p1.globalIndex), incUR.get(p2.globalIndex),
//				    					incUR.get(p3.globalIndex));
//			    			}
			    			rlt.set(i, r);
			    		}
		    		}
	    		}
	    	}
	    }
	    return rlt;	
	}
	
	public void run(int elementType, String gridFile, String outputFolder) {
		MeshReader reader = null;
		reader = new MeshReader(gridFile);
		Mesh mesh = reader.read2DMesh();
		
		this.outputFolder = outputFolder;
		
		ScalarShapeFunction[] shapeFun = null;
		ScalarShapeFunction[] shapeFunRect = null;
		if(elementType == 1) {
			mesh.computeNodesBelongToElement();
			mesh.computeNeiborNode();
			//Assign degree of freedom to element
			shapeFun = new SFLinearLocal2D[3];
			for(int i=0;i<3;i++)
				shapeFun[i] = new SFLinearLocal2D(i+1);
			
			shapeFunRect = new SFBilinearLocal2D[4];
			for(int i=0;i<4;i++)
				shapeFunRect[i] = new SFBilinearLocal2D(i+1);
			
			//Assign shape function to DOF
			for(int i=1;i<=mesh.getElementList().size();i++) {
				Element e = mesh.getElementList().at(i);
				if(e.nodes.size() == 4) {
					for(int j=1;j<=e.nodes.size();j++) {
						DOF dof = new DOF(j,e.nodes.at(j).globalIndex,shapeFunRect[j-1]);
						e.addNodeDOF(j, dof);
					}
				} else if(e.nodes.size() == 3) {
					for(int j=1;j<=e.nodes.size();j++) {
						DOF dof = new DOF(j,e.nodes.at(j).globalIndex,shapeFun[j-1]);
						e.addNodeDOF(j, dof);
					}
				} else {
					System.out.println("Error: e.nodes.size()="+e.nodes.size());
				}
			}
		} else if(elementType == 2) {
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
				
//				for(int j=1;j<=3;j++) {
//					Node l = e.nodes.at(indexSet[j-1]);
//					Node r = e.nodes.at(indexSet[j]);
//					double cx = (l.coord(1)+r.coord(1))/2.0;
//					double cy = (l.coord(2)+r.coord(2))/2.0;
//					Node node = new Node(2);
//					node.set(mesh.getNodeList().size()+1, cx,cy );
//					Node findNode = mesh.containNode(node);
//					if(findNode == null) {
//						e.addVertexNode(node, false);
//						mesh.addNode(node);
//					} else {
//						e.addVertexNode(findNode, false);
//					}
//				}
				
				
			}
			mesh.computeNodesBelongToElement();
			mesh.computeNeiborNode();
			//Assign degree of freedom to element
//			shapeFun = new SFQuadraticLocal2DFast[6];
//			for(int i=1;i<=6;i++)
//				shapeFun[i-1] = new SFQuadraticLocal2DFast(i);
			shapeFun = new SFQuadraticLocal2D[6];
			for(int i=1;i<=6;i++)
				shapeFun[i-1] = new SFQuadraticLocal2D(i);

			//Assign shape function to DOF
			for(int i=1;i<=mesh.getElementList().size();i++) {
				Element e = mesh.getElementList().at(i);
				for(int j=1;j<=e.nodes.size();j++) {
					DOF dof = new DOF(j,e.nodes.at(j).globalIndex,shapeFun[j-1]);
					e.addNodeDOF(j, dof);
				}
			}
		} else {
			System.out.println("Error: elementType parameter!");
			return;
		}


		
		NodeList list = mesh.getNodeList();
		int nNode = list.size();
		
		//----------------------Begin Left light source ---------------------
		setDelta(1.0, 2.8);
		
		//Solve background forward problem
		setMu_a(0.0, 0.0, 0.0, 
				0.1, 1);
		Vector bkUL = solveForwardNeumann(mesh);
		plotVector(mesh, bkUL, "prostate_test1_bkUL.dat");
		
		//Solve forward problem with inclusion
		setMu_a(2.0, 2.6, 0.3, 
				1.0, 1);
		plotFunction(mesh, this.mu_a, "prostate_test1_alpha_real.dat");
		Vector incUL = solveForwardNeumann(mesh);
		plotVector(mesh, incUL, "prostate_test1_incUL.dat");
		//中间不要加入其他代码！
	    Vector tailUL = computeTailLeftLightSource(mesh, bkUL, incUL);
	    plotVector(mesh, tailUL, "prostate_test1_tailUL.dat");
	    
		//Cut the peak of incUL, so that we can see the signification counter line at the position of inclusion
	    Vector tailUL_cut = new SparseVector(tailUL.getDim());
	    for(int i=1;i<=tailUL.getDim();i++)
	    	if(tailUL.get(i)<80000)
	    		tailUL_cut.set(i, tailUL.get(i));
	    plotVector(mesh, tailUL_cut, "prostate_test1_tailUL_cut.dat");
//	    tailUL = Utils.gaussSmooth(mesh, tailUL, 1, 0.5);
//	    tailUL = Utils.gaussSmooth(mesh, tailUL, 1, 0.5);
//	    tailUL = Utils.gaussSmooth(mesh, tailUL, 1, 0.5);
//	    tailUL = Utils.gaussSmooth(mesh, tailUL, 1, 0.5);
//	    plotVector(mesh, tailUL, "prostate_test1_tailUL_smooth.dat");		
	    
	    //incU_x
		plotVector(mesh, Tools.computeDerivative(mesh,incUL,"x"), "prostate_test1_incUL_x.dat");
	    
		//Recovery parameter mu_a from solution
	    Vector alpha_real_cmp = solveParamInverse(mesh,incUL);
	    plotVector(mesh, alpha_real_cmp, "prostate_test1_alpha_real_cmp.dat");

	    //alpha_real_cmp_smooth_x
	    Vector alpha_real_cmp_smooth = Utils.gaussSmooth(mesh, alpha_real_cmp, 1, 0.5);
	    alpha_real_cmp_smooth = Utils.gaussSmooth(mesh, alpha_real_cmp_smooth, 1, 0.5);
		plotVector(mesh, Tools.computeDerivative(mesh,alpha_real_cmp_smooth,"x"), "prostate_test1_alpha_real_cmp_x.dat");

		//Cut the peak of incUL, so that we can see the signification counter line at the position of inclusion
	    Vector incUL_cut = new SparseVector(incUL.getDim());
	    for(int i=1;i<=incUL.getDim();i++)
	    	if(incUL.get(i)<80000.0)
	    		incUL_cut.set(i, incUL.get(i));
	    plotVector(mesh, incUL_cut, "prostate_test1_incUL_cut.dat");
	    
	    
	    //Recovery parameter mu_a from tail
	    Vector tailUL_noise = addNoise(tailUL,0.0);
	    Vector alphaL = solveParamInverse(mesh,tailUL_noise);
	    //Cut noise
	    for(int i=1;i<=nNode;i++) {
	    	Node node = list.at(i);
	    	if(node.coord(1) <1.3 || node.coord(1)>4.4 || node.coord(2)<2.0) {
		    //if(node.coord(1) <1.9 || node.coord(1)>2.8 || node.coord(2)<1.5) {
	    		alphaL.set(node.globalIndex, 0.1);
	    	}
		}
	    plotVector(mesh, alphaL, "prostate_test1_alphaL.dat");
	    	    
		//----------------------Begin Right light source ---------------------
		setDelta(4.0, 2.8);
		
		setMu_a(0.0, 0.0, 0.0, 
				0.1, 1);
		Vector bkUR = solveForwardNeumann(mesh);
		plotVector(mesh, bkUR, "prostate_test1_bkUR.dat");
		
		setMu_a(2.0, 2.6, 0.3, 
				1.0, 1);
		Vector incUR = solveForwardNeumann(mesh);
	    plotVector(mesh, incUR, "prostate_test1_incUR.dat");
	    Vector tailUR = computeTailRightLightSource(mesh, bkUR, incUR);
	    plotVector(mesh, tailUR, "prostate_test1_tailUR.dat");
//	    tailUR = Utils.gaussSmooth(mesh, tailUR, 1, 0.5);
//	    tailUR = Utils.gaussSmooth(mesh, tailUR, 1, 0.5);
//	    tailUR = Utils.gaussSmooth(mesh, tailUR, 1, 0.5);
//	    tailUR = Utils.gaussSmooth(mesh, tailUR, 1, 0.5);
//	    plotVector(mesh, tailUR, "prostate_test1_tailUR_smooth.dat");
	    
	    Vector alphaR = solveParamInverse(mesh,tailUR);
	    //Cut noise
	    for(int i=1;i<=nNode;i++) {
	    	Node node = list.at(i);
	    	if(node.coord(1) <0.6 || node.coord(1)>3.7 || node.coord(2)<2.0) {
		    //if(node.coord(1) <1.9 || node.coord(1)>2.8 || node.coord(2)<1.5) {
	    		alphaR.set(node.globalIndex, 0.1);
	    	}
		}
	    plotVector(mesh, alphaR, "prostate_test1_alphaR.dat");
	    
	    //Smooth...
	    alphaL = Utils.gaussSmooth(mesh, alphaL, 1, 0.5);
	    alphaL = Utils.gaussSmooth(mesh, alphaL, 1, 0.5);
	    alphaL = Utils.gaussSmooth(mesh, alphaL, 1, 0.5);
	    alphaR = Utils.gaussSmooth(mesh, alphaR, 1, 0.5);
	    alphaR = Utils.gaussSmooth(mesh, alphaR, 1, 0.5);
	    alphaR = Utils.gaussSmooth(mesh, alphaR, 1, 0.5);

	    Vector alpha_avg = new SparseVector(nNode);
	    for(int i=1;i<=nNode;i++) {
	    	alpha_avg.set(i, (alphaL.get(i)+alphaR.get(i))/2.0);
	    }
	    plotVector(mesh, alpha_avg, "prostate_test1_alpha_avg.dat");
	    
	    //Smooth...
	    alpha_avg = Utils.gaussSmooth(mesh, alpha_avg, 2, 0.5);
	    alpha_avg = Utils.gaussSmooth(mesh, alpha_avg, 2, 0.5);
	    alpha_avg = Utils.gaussSmooth(mesh, alpha_avg, 2, 0.5);
	    Vector alpha_avg_smooth = Utils.gaussSmooth(mesh, alpha_avg, 2, 0.5);
	    plotVector(mesh, alpha_avg_smooth, "prostate_test1_alpha_avg_smooth.dat");
	    
	    Double max = alpha_avg_smooth.normInf();
	    for(int i=1;i<=alpha_avg_smooth.getDim();i++) {
	    	if(alpha_avg_smooth.get(i) < 0.8*max)
	    		alpha_avg_smooth.set(i, 0.1);
	    }
	    plotVector(mesh, alpha_avg_smooth, "prostate_test1_alpha_avg_smooth_no_noise.dat");
	    
	    
	    //-------------------- Begin Nonlinear Enhancement -----------------------

	    Vector alpha_m1 = alpha_avg_smooth;
		setDelta(1.0, 2.8);
		setMu_a(2.0, 2.6, 0.3, 
				1.0, 1);
		//！！！！！！！！！mu_a = new VectorBasedFunction(alpha_m1);
		//???边界条件需要满足测量值，而不是Neumann边界
		Function diri = new Vector2Function(incUL);
		Vector um1 = solveForwardDirichlet(mesh,diri);
		//Vector um1 = solveForwardNeumann(mesh);
		plotVector(mesh, um1, "prostate_test1_enhance_um1.dat");
		
//		//Nonlinear iteration
//		for(int nit=2;nit<=6;nit++) {
//		    Vector alpha_m2 = solveParamInverse(mesh,um1);
//			plotVector(mesh, alpha_m2, "prostate_test1_enhance_alpha_m"+nit+".dat");
//	
//		    Vector wL1 = solveEnhance(mesh,alpha_m2,alpha_m1,um1,nit);
//		    plotVector(mesh, wL1, "prostate_test1_enhance_wL"+(nit-1)+".dat");
//		    
//		    Vector um2 = Vector.axpy(1.0, um1, wL1);
//		    plotVector(mesh, um2, "prostate_test1_enhance_um"+nit+".dat");
//		    
//		    um1 = um2;
//		    alpha_m1 = alpha_m2;
//		}
		
		//Global Convergence Method
		int N=3;
		double[]s = new double[N];
		double h = 0.1;
		s[0] = 2.0;
		// 由远到近，s递减
		for(int i=1; i<N; i++)
			s[i] = s[0] - i*h;
		Vector tailT = new SparseVector(um1.getDim());
		// tail函数（最远处）
		for(int i=1;i<=tailT.getDim();i++) {
			tailT.set(i, Math.log(um1.get(i))/(s[0]*s[0]));
		}
		this.plotVector(mesh, tailT, "tailT.dat");
		
		//Simulate border measurement data
		Vector[] ui = new Vector[N];
		Vector[] phi = new Vector[N];
		setMu_a(2.0, 2.6, 0.3, 
				1.0, 1);
		//mu_a = new VectorBasedFunction(alpha_m1);
		for(int i=0;i<N;i++) {
			setDelta(1.0+i*h, 2.8);
			ui[i] = solveForwardNeumann(mesh);
			plotVector(mesh, ui[i], "GCMincU_"+i+".dat");
		}
		for(int i=0;i<N-1;i++) {
			int dim = ui[i].getDim();
			phi[i] = new SparseVector(dim);
			for(int j=1;j<=dim;j++) {
				phi[i].set(j, 
						(1.0/h) * 
						(Math.log(ui[i  ].get(j))/(s[i  ]*s[i  ]) - 
						 Math.log(ui[i+1].get(j))/(s[i+1]*s[i+1]) 
						));
			}
			plotVector(mesh, phi[i], "GCMphi_"+i+".dat");
		}
		
//		//模拟计算出来的phi除了提供边界条件，区域内部的值也有，可以用来重构真实的a(x)
//		int dim=um1.getDim();
//		Vector v_tidle_real = new Vector(dim);
//		for(int j=0;j<N-1;j++) {
//			v_tidle_real = Vector.axpy(-s[N-1]*s[N-1]*h, phi[j], v_tidle_real);
//		}
//		v_tidle_real = Vector.axpy(s[N-1]*s[N-1], tailT, v_tidle_real);
//		plotVector(mesh, v_tidle_real, "v_tidle_N_real.dat");
//		
//		Vector u_real = new Vector(dim);
//		for(int i=1;i<=dim;i++) {
//			u_real.set(i, Math.pow(Math.E, v_tidle_real.get(i)));
//		}
//		plotVector(mesh, u_real, "u_N_real.dat");
//		
//		Vector a_real = solveParamInverse(mesh, u_real);
//		plotVector(mesh, a_real, "alpha_recon_N_real.dat");
		
		
		solveGCM(mesh, N, s, phi, tailT);
		
	}
	
	public Vector addNoise(Vector v, double persent) {
		Vector rlt = new SparseVector(v.getDim());
		for(int i=1;i<=v.getDim();i++) {
			double val = v.get(i);
			val += val*persent*(2*Math.random()-1.0);
			rlt.set(i, val);
		}
		return rlt;
	}
	
	public void assignLinearShapFunction(Mesh mesh) {
		ScalarShapeFunction[] shapeFun = null;
		ScalarShapeFunction[] shapeFunHalf = null;
		ScalarShapeFunction[] shapeFunRect = null;
		ScalarShapeFunction[] shapeFunRectHalf = null;
		
		//Assign degree of freedom to element
		shapeFun = new SFLinearLocal2D[3];
		shapeFunHalf = new SFLinearLocal2D[3];
		for(int i=0;i<3;i++) {
			shapeFun[i] = new SFLinearLocal2D(i+1);
			shapeFunHalf[i] = new SFLinearLocal2D(i+1,0.5);
		}
		
		shapeFunRect = new SFBilinearLocal2D[4];
		shapeFunRectHalf = new SFBilinearLocal2D[4];
		for(int i=0;i<4;i++) {
			shapeFunRect[i] = new SFBilinearLocal2D(i+1);
			shapeFunRectHalf[i] = new SFBilinearLocal2D(i+1,0.5);
		}
		
		//Assign shape function to DOF
		for(int i=1;i<=mesh.getElementList().size();i++) {
			Element e = mesh.getElementList().at(i);
			e.clearAllDOF();
			if(e.nodes.size() == 4) {
				//Asign degree of freedom to element
				int nDofLocalIndexCounter = 0;
				for(int j=1;j<=e.nodes.size();j++) {
					//Asign shape function to DOF
					if(e.nodes.at(j) instanceof NodeRefined) {
						NodeRefined nRefined = (NodeRefined)e.nodes.at(j);
						if(nRefined.isHangingNode()) {
							DOF dof  = new DOF(++nDofLocalIndexCounter,nRefined.constrainNodes.at(1).globalIndex,
									shapeFunRectHalf[j-1]);
							e.addNodeDOF(j, dof);
							DOF dof2 = new DOF(++nDofLocalIndexCounter,nRefined.constrainNodes.at(2).globalIndex,
									shapeFunRectHalf[j-1]);
							e.addNodeDOF(j, dof2);
						} else {
							DOF dof = new DOF(++nDofLocalIndexCounter,e.nodes.at(j).globalIndex,shapeFunRect[j-1]);
							e.addNodeDOF(j, dof);				
						}
					} else {
						DOF dof = new DOF(++nDofLocalIndexCounter,e.nodes.at(j).globalIndex,shapeFunRect[j-1]);
						e.addNodeDOF(j, dof);
					}
				}
			} else if(e.nodes.size() == 3) {
				int nDofLocalIndexCounter = 0;
				for(int j=1;j<=e.nodes.size();j++) {
					//Asign shape function to DOF
					if(e.nodes.at(j) instanceof NodeRefined) {
						NodeRefined nRefined = (NodeRefined)e.nodes.at(j);
						if(nRefined.isHangingNode()) {
							DOF dof  = new DOF(++nDofLocalIndexCounter,nRefined.constrainNodes.at(1).globalIndex,
									shapeFunHalf[j-1]);
							e.addNodeDOF(j, dof);
							DOF dof2 = new DOF(++nDofLocalIndexCounter,nRefined.constrainNodes.at(2).globalIndex,
									shapeFunHalf[j-1]);
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
				
			} else {
				System.out.println("Error: e.nodes.size()="+e.nodes.size());
			}
		}
	}
	
	
	
	public Vector solveAdaptive(Mesh mesh, int elementType, String outputFolder) {
		this.outputFolder = outputFolder;

		if(elementType == 1) {
			mesh.computeNodesBelongToElement();
			mesh.computeNeiborNode();
			mesh.computeNeighborElement();
			assignLinearShapFunction(mesh);
		} else {
			System.out.println("Error: elementType parameter!");
			return null;
		}

		NodeList list = mesh.getNodeList();
		int nNode = list.size();
		
		//----------------------Begin Left light source ---------------------
		setDelta(1.0, 2.8);
		
		//Solve background forward problem
		setMu_a(0.0, 0.0, 0.0, 
				0.1, 1);
		Vector bkUL = solveForwardNeumann(mesh);
		plotVector(mesh, bkUL, "prostate_test1_bkUL.dat");
		
		//  \Laplace U
		//plotVector(mesh, computeLaplace2D(mesh, bkUL), "prostate_test1_bkUL_laplace.dat");
		
		//Solve forward problem with inclusion
		setMu_a(2.5, 2.6, 0.3, 
				0.2, 1);
		plotFunction(mesh, this.mu_a, "prostate_test1_alpha_real.dat");
		Vector incUL = solveForwardNeumann(mesh);
		plotVector(mesh, incUL, "prostate_test1_incUL.dat");
		//中间不要加入其他代码！
	    Vector tailUL = computeTailLeftLightSource(mesh, bkUL, incUL);
	    plotVector(mesh, tailUL, "prostate_test1_tailUL.dat");
	    plotTopBorder(mesh,incUL, "prostate_test1_incUL_topLine.dat");
	    
	    
	    
	    //incU_x
		//plotVector(mesh, computeDerivative(mesh,incUL,"x"), "prostate_test1_incUL_x.dat");
		//Recovery parameter mu_a from solution
	    Vector alpha_real_cmp = solveParamInverse(mesh,incUL);
	    plotVector(mesh, alpha_real_cmp, "prostate_test1_alpha_real_cmp.dat");

	    //alpha_real_cmp_smooth_x
	    Vector alpha_real_cmp_x = Tools.computeDerivative(mesh,alpha_real_cmp,"x");
		plotVector(mesh, alpha_real_cmp_x, "prostate_test1_alpha_real_cmp_x.dat");

	    //Recovery parameter mu_a from tail
	    Vector tailUL_noise = addNoise(tailUL,0.0);
	    Vector alphaL = solveParamInverse(mesh,tailUL_noise);
	    plotVector(mesh, alphaL, "prostate_test1_alphaL_before_cut.dat");
	    //Cut noise
	    for(int i=1;i<=nNode;i++) {
	    	Node node = list.at(i);
	    	//if(node.coord(1) <1.3 || node.coord(1)>4.4 || node.coord(2)<2.0) {
		    //if(node.coord(1) <1.9 || node.coord(1)>2.8 || node.coord(2)<1.5) {
		    //incusion: (2.5,2.75),r=0.15
			if(node.coord(1) <1.5 || node.coord(1)>3.5 || node.coord(2)<2.0) {
	    		alphaL.set(node.globalIndex, 0.1);
	    	}
		}
	    plotVector(mesh, alphaL, "prostate_test1_alphaL.dat");
	    alphaL = Utils.gaussSmooth(mesh, alphaL, 1, 0.8);
	    alphaL = Utils.gaussSmooth(mesh, alphaL, 1, 0.8);
	    alphaL = Utils.gaussSmooth(mesh, alphaL, 1, 0.8);
	    alphaL = Utils.gaussSmooth(mesh, alphaL, 1, 0.8);
	    plotVector(mesh, alphaL, "prostate_test1_alphaL_smooth.dat");

	    	    
		//----------------------Begin Right light source ---------------------
		setDelta(4.0, 2.8);
		
		setMu_a(0.0, 0.0, 0.0, 
				0.1, 1);
		Vector bkUR = solveForwardNeumann(mesh);
		plotVector(mesh, bkUR, "prostate_test1_bkUR.dat");
		
		setMu_a(2.5, 2.75, 0.15, 
				2.0, 1);
		Vector incUR = solveForwardNeumann(mesh);
	    plotVector(mesh, incUR, "prostate_test1_incUR.dat");
	    Vector tailUR = computeTailRightLightSource(mesh, bkUR, incUR);
	    plotVector(mesh, tailUR, "prostate_test1_tailUR.dat");
	    
	    Vector alphaR = solveParamInverse(mesh,tailUR);
	    plotVector(mesh, alphaR, "prostate_test1_alphaR_before_cut.dat");
	    //Cut noise
	    for(int i=1;i<=nNode;i++) {
	    	Node node = list.at(i);
	    	//if(node.coord(1) <0.6 || node.coord(1)>3.7 || node.coord(2)<2.0) {
		    //if(node.coord(1) <1.9 || node.coord(1)>2.8 || node.coord(2)<1.5) {
	    	//incusion: (2.5,2.75),r=0.15
			if(node.coord(1) <1.5 || node.coord(1)>3.5 || node.coord(2)<2.0) {
	    		alphaR.set(node.globalIndex, 0.1);
	    	}
		}
	    plotVector(mesh, alphaR, "prostate_test1_alphaR.dat");
	    alphaR = Utils.gaussSmooth(mesh, alphaR, 1, 0.8);
	    alphaR = Utils.gaussSmooth(mesh, alphaR, 1, 0.8);
	    alphaR = Utils.gaussSmooth(mesh, alphaR, 1, 0.8);
	    alphaR = Utils.gaussSmooth(mesh, alphaR, 1, 0.8);
	    plotVector(mesh, alphaR, "prostate_test1_alphaR_smooth.dat");
	    
	    Vector alpha_avg = new SparseVector(nNode);
	    for(int i=1;i<=nNode;i++) {
	    	alpha_avg.set(i, (alphaL.get(i)+alphaR.get(i))/2.0);
	    }
	    plotVector(mesh, alpha_avg, "prostate_test1_alpha_avg.dat");
	    
	    return alpha_avg;
	}
	
	
	
	public void runAdaptive(int elementType, int gridID, String outputFolder) {
		MeshReader reader = null;
		reader = new MeshReader("prostate_test"+gridID+".grd");
		Mesh mesh = reader.read2DMesh();

		Vector alpha_avg = solveAdaptive(mesh, elementType, outputFolder);
	    //Smooth...
		Vector alpha_avg_smooth = Utils.gaussSmooth(mesh, alpha_avg, 1, 0.5);
	    alpha_avg_smooth = Utils.gaussSmooth(mesh, alpha_avg_smooth, 1, 0.5);
	    alpha_avg_smooth = Utils.gaussSmooth(mesh, alpha_avg_smooth, 1, 0.5);
	    alpha_avg_smooth = Utils.gaussSmooth(mesh, alpha_avg_smooth, 1, 0.5);
	    alpha_avg_smooth = Utils.gaussSmooth(mesh, alpha_avg_smooth, 1, 0.5);
	    alpha_avg_smooth = Utils.gaussSmooth(mesh, alpha_avg_smooth, 1, 0.5);
	    plotVector(mesh, alpha_avg_smooth, "prostate_test1_alpha_avg_smooth.dat");
		Vector alpha_avg_smooth_no_noise = alpha_avg_smooth.copy();
	    Double max = alpha_avg_smooth_no_noise.normInf();
	    for(int i=1;i<=alpha_avg_smooth_no_noise.getDim();i++) {
	    	if(alpha_avg_smooth_no_noise.get(i) < 0.7*max)
	    		alpha_avg_smooth_no_noise.set(i, 0.1);
	    }
	    plotVector(mesh, alpha_avg_smooth_no_noise, "prostate_test1_alpha_avg_smooth_no_noise.dat");		
		
		//Adaptive refinement 1
		ElementList eToRefine = computeRefineElement(mesh, alpha_avg_smooth, 0.06);
		System.out.println("Before refine: Element="+mesh.getElementList().size()+", Node="+mesh.getNodeList().size());
		Refiner.refineOnce(mesh, eToRefine);
		System.out.println("After refine: Element="+mesh.getElementList().size()+", Node="+mesh.getNodeList().size());
		plotFunction(mesh, this.mu_a, "prostate_test1_alpha_real_refine.dat");
	
		alpha_avg = solveAdaptive(mesh, elementType, outputFolder+"_adaptive");
		alpha_avg_smooth = Utils.gaussSmooth(mesh, alpha_avg, 1, 0.5);
		alpha_avg_smooth = Utils.gaussSmooth(mesh, alpha_avg, 1, 0.5);
		alpha_avg_smooth = Utils.gaussSmooth(mesh, alpha_avg, 1, 0.5);
		alpha_avg_smooth = Utils.gaussSmooth(mesh, alpha_avg, 1, 0.5);
	    alpha_avg_smooth = Utils.gaussSmooth(mesh, alpha_avg, 1, 0.5);
	    alpha_avg_smooth = Utils.gaussSmooth(mesh, alpha_avg, 1, 0.5);
	    plotVector(mesh, alpha_avg_smooth, "prostate_test1_alpha_avg_smooth.dat");
	    alpha_avg_smooth_no_noise = alpha_avg_smooth.copy();
	    max = alpha_avg_smooth_no_noise.normInf();
	    for(int i=1;i<=alpha_avg_smooth_no_noise.getDim();i++) {
	    	if(alpha_avg_smooth_no_noise.get(i) < 0.7*max)
	    		alpha_avg_smooth_no_noise.set(i, 0.1);
	    }
	    plotVector(mesh, alpha_avg_smooth_no_noise, "prostate_test1_alpha_avg_smooth_no_noise.dat");		
		
		//Adaptive refinement 2
		eToRefine.clear();
		eToRefine = computeRefineElement(mesh, alpha_avg_smooth, 0.10);
		System.out.println("Before refine: Element="+mesh.getElementList().size()+", Node="+mesh.getNodeList().size());
		Refiner.refineOnce(mesh, eToRefine);
		System.out.println("After refine: Element="+mesh.getElementList().size()+", Node="+	mesh.getNodeList().size());
		plotFunction(mesh, this.mu_a, "prostate_test1_alpha_real_refine.dat");
		
		alpha_avg = solveAdaptive(mesh, elementType, outputFolder+"_adaptive2");
	    alpha_avg_smooth = Utils.gaussSmooth(mesh, alpha_avg, 1, 0.5);
	    alpha_avg_smooth = Utils.gaussSmooth(mesh, alpha_avg, 1, 0.5);
	    alpha_avg_smooth = Utils.gaussSmooth(mesh, alpha_avg, 1, 0.5);
	    plotVector(mesh, alpha_avg_smooth, "prostate_test1_alpha_avg_smooth.dat");
	    alpha_avg_smooth_no_noise = alpha_avg_smooth.copy();
	    max = alpha_avg_smooth_no_noise.normInf();
	    for(int i=1;i<=alpha_avg_smooth_no_noise.getDim();i++) {
	    	if(alpha_avg_smooth_no_noise.get(i) < 0.7*max)
	    		alpha_avg_smooth_no_noise.set(i, 0.1);
	    }
	    plotVector(mesh, alpha_avg_smooth_no_noise, "prostate_test1_alpha_avg_smooth_no_noise.dat");		
	}
	
	public ElementList computeRefineElement(Mesh mesh, Vector v, double persent) {
	    Vector v_smooth = Utils.gaussSmooth(mesh, v, 1, 0.5);

		ElementList eList = mesh.getElementList();
		ElementList eToRefine = new ElementList();

		List<PairDoubleInteger> list = new ArrayList<PairDoubleInteger>();
		for(int i=1;i<=eList.size();i++) {
			Element e = eList.at(i);
			Vector v1 = new SparseVector(e.nodes.size());
			for(int j=1;j<=e.nodes.size();j++) {
				Node node = e.nodes.at(j);
				v1.set(j, v.get(node.globalIndex)-v_smooth.get(node.globalIndex));
			}
			list.add(new PairDoubleInteger(v1.norm2(),i));
		}
		Collections.sort(list, new Comparator<PairDoubleInteger>() {
			@Override
			public int compare(PairDoubleInteger o1, PairDoubleInteger o2) {
				return o1.d < o2.d ? 1 : -1;
			}
		});
		int nThreshold = (int) Math.floor(persent*eList.size());
		for(int i=1;i<=nThreshold;i++) {
			eToRefine.add(eList.at(list.get(i).i));
		}
		return eToRefine;
	}

	public void plotTopBorder(Mesh mesh, Vector v, String fileName) {
		NodeList nList = mesh.getNodeList();
		int nNode = nList.size();
		List<PairDoubleInteger> pList = new ArrayList<PairDoubleInteger>();
		for(int i=1;i<=nNode;i++) {
			if(Math.abs(nList.at(i).coord(2)-3.0)<Constant.eps) {
				PairDoubleInteger pair = new PairDoubleInteger(nList.at(i).coord(1), i);
				pList.add(pair);
			}
		}
		Collections.sort(pList, new Comparator<PairDoubleInteger>(){
			@Override
			public int compare(PairDoubleInteger o1, PairDoubleInteger o2) {
				return o1.d>o2.d?1:-1;
			}
		});
	    MeshWriter writer = new MeshWriter(mesh);
	    if(!outputFolder.isEmpty()) {
		    File file = new File("./"+outputFolder);
			if(!file.exists()) {
				file.mkdirs();
			}
	    }
	    Vector x = new SparseVector(pList.size());
	    Vector y = new SparseVector(pList.size());
	    for(int i=0;i<pList.size();i++) {
	    	PairDoubleInteger pair = pList.get(i);
	    	x.set(i+1, pair.d);
	    	y.set(i+1, v.get(pair.i));
	    }
	    writer.writeTechplotLine("./"+outputFolder+"/"+fileName, x,y);
		
		
	}
	
	/**
	 * 
	 * @param mesh1
	 * @param mesh2
	 * @param u
	 * @return
	 */
	public Vector extractData(Mesh meshFrom, Mesh meshTo, Vector u) {
		NodeList nodeTo = meshTo.getNodeList();
		int dimTo = nodeTo.size();
		Vector rlt = new SparseVector(dimTo);
		for(int i=1;i<=nodeTo.size();i++) {
			Node node = meshFrom.containNode(nodeTo.at(i));
			if(node != null) {
				rlt.set(i, u.get(node.globalIndex));
			}
		}
		return rlt;
	}	
	
	Vector getGCMCoef(double sn_1, double sn) {
		Vector A = new SparseVector(4);
		double I0 =  3*(sn_1-sn) - 2*sn_1*(Math.log(sn_1)-Math.log(sn));
		double A1 = (1.0/I0)*(
				(-3.0/2.0)*(Math.pow(sn_1, 4)-Math.pow(sn, 4))+
				(10.0/3.0)*sn_1*(Math.pow(sn_1, 3)-Math.pow(sn, 3))+
				(-2.0)*sn_1*sn_1*(Math.pow(sn_1, 2)-Math.pow(sn, 2))
				);
		double A2 = (1.0/I0)*(
				(-10.0/3.0)*(Math.pow(sn_1, 3)-Math.pow(sn, 3))+
				4.0*sn_1*(Math.pow(sn_1, 2)-Math.pow(sn, 2))
				);
		double A3 = (-2.0/I0)*(Math.log(sn_1)-Math.log(sn));
		double A4 = (-2.0/I0)*(Math.pow(sn_1, 2)-Math.pow(sn, 2));
		A.set(1, A1);
		A.set(2, A2);
		A.set(3, A3);
		A.set(4, A4);
		return A;
	}
	
	/**
	 * Test for Global Convergence Method 
	 * 下一步需要做的事情：
	 * 光源s从最远到最近变化，最远距离选一个足够远的距离，h足够小，迭代足够多的步数后，
	 * 观察Tail的影响有多大？是否想之前的测试结果，影响很大？
	 * 
	 * @param elementType
	 * @param gridFileForward
	 * @param gridFileGCM
	 * @param outputFolder
	 */
	public void runGCMTest(int elementType, 
			String gridFileForward, 
			String gridFileGCM, 
			String outputFolder) {
		MeshReader readerForward = new MeshReader(gridFileForward);
		Mesh meshForward = readerForward.read2DMesh();
		MeshReader readerGCM = new MeshReader(gridFileGCM);
		Mesh meshGCM = readerGCM.read2DMesh();
		
		this.outputFolder = outputFolder;
		
		ScalarShapeFunction[] shapeFun = null;
		ScalarShapeFunction[] shapeFunRect = null;
		if(elementType == 1) {
			//Assign degree of freedom to element
			shapeFun = new SFLinearLocal2D[3];
			for(int i=0;i<3;i++)
				shapeFun[i] = new SFLinearLocal2D(i+1);
			
			shapeFunRect = new SFBilinearLocal2D[4];
			for(int i=0;i<4;i++)
				shapeFunRect[i] = new SFBilinearLocal2D(i+1);
			
			meshForward.computeNodesBelongToElement();
			meshForward.computeNeiborNode();
			//Assign shape function to DOF
			for(int i=1;i<=meshForward.getElementList().size();i++) {
				Element e = meshForward.getElementList().at(i);
				if(e.nodes.size() == 4) {
					for(int j=1;j<=e.nodes.size();j++) {
						DOF dof = new DOF(j,e.nodes.at(j).globalIndex,shapeFunRect[j-1]);
						e.addNodeDOF(j, dof);
					}
				} else if(e.nodes.size() == 3) {
					for(int j=1;j<=e.nodes.size();j++) {
						DOF dof = new DOF(j,e.nodes.at(j).globalIndex,shapeFun[j-1]);
						e.addNodeDOF(j, dof);
					}
				} else {
					System.out.println("Error: e.nodes.size()="+e.nodes.size());
				}
			}
			meshGCM.computeNodesBelongToElement();
			meshGCM.computeNeiborNode();
			//Assign shape function to DOF
			for(int i=1;i<=meshGCM.getElementList().size();i++) {
				Element e = meshGCM.getElementList().at(i);
				if(e.nodes.size() == 4) {
					for(int j=1;j<=e.nodes.size();j++) {
						DOF dof = new DOF(j,e.nodes.at(j).globalIndex,shapeFunRect[j-1]);
						e.addNodeDOF(j, dof);
					}
				} else if(e.nodes.size() == 3) {
					for(int j=1;j<=e.nodes.size();j++) {
						DOF dof = new DOF(j,e.nodes.at(j).globalIndex,shapeFun[j-1]);
						e.addNodeDOF(j, dof);
					}
				} else {
					System.out.println("Error: e.nodes.size()="+e.nodes.size());
				}
			}			

		} else {
			System.out.println("Error: elementType parameter!");
			return;
		}

		NodeList list = meshForward.getNodeList();
		int nNode = list.size();

		double coefK = 0.02;
		
		setDelta(1.0, 3.5);
		//setDelta(0.0, 4.0);
		plotFunction(meshForward, this.delta, "delta.dat");
		
		//Solve background forward problem
		setMu_a(0.0, 0.0, 0.0, 
				0.1, 1);
		Vector bkUL = solveForwardNeumann(meshForward);
		plotVector(meshForward, bkUL, "bkUL.dat");
		
		//Solve forward problem with inclusion
		setMu_a(2.0, 2.9, 0.3, 
				2.0, 1);
		plotFunction(meshForward, this.mu_a, "alpha_real.dat");
		Vector incUL = solveForwardNeumann(meshForward);
		plotVector(meshForward, incUL, "incUL.dat");
		
		Vector incUL_x = Tools.computeDerivative(meshForward, incUL, "x");
		plotVector(meshForward, incUL_x, "incUL_x.dat");
		Vector incUL_xx = Tools.computeDerivative(meshForward, incUL_x, "x");
		plotVector(meshForward, incUL_xx, "incUL_xx.dat");
		
		Vector incUL_laplace = Tools.computeLaplace2D(meshForward, incUL);
		plotVector(meshForward, incUL_laplace, "incUL_laplace.dat");
		
		Vector vMu_a = discreteFunction(meshForward,this.mu_a);
		Vector eq = SparseVector.axMuly(1.0, vMu_a, incUL);
		plotVector(meshForward, eq, "incUL_Mu_amU.dat");
		eq = SparseVector.axpy(-coefK, incUL_laplace, eq);
		plotVector(meshForward, eq, "incUL_equation.dat");


		//---------------Global Convergence Method---------------
		//Number of moving light sources
		int N=4; 
		double[]s = new double[N];
		double h = 0.4; 
		//h越小越好？（由于错误的代码(sn-1 - sn)，h=0.01结果很好，再大结果就很坏）
		//正确代码：h=0.1没问题,h=0.05要好一些
		s[0] = 3.0; //s[0]>=1?  =1 NO; =2,4 OK
		for(int i=1; i<N; i++)
			s[i] = s[0] - i*h; //由远到近，s递减
		
		//Simulated border measurement data
		int dimGCM = meshGCM.getNodeList().size();
		Vector[] ui_ex = new Vector[N];
		Vector[] ui = new Vector[N];
		Vector[] vi = new Vector[N];
		Vector[] phi = new Vector[N];
		
		//ui[0] ui[1] ui[2]
		//vi[0] vi[1] vi[2]
		for(int i=0;i<N;i++) {
			//实际上，应该来自从tail重构出来的a(x)
			setMu_a(2.0, 2.9, 0.3, 
					2.0, 1);
			
			setDelta(1.0+i*h, 3.5);
			
			ui_ex[i] = solveForwardNeumann(meshForward); //在meshGCM上求解，top边界条件不好确定
			plotVector(meshForward, ui_ex[i], "GCM_incU_ex_"+i+".dat");
			
			//截取meshForward的部分解到meshGCM上
			ui[i] = extractData(meshForward, meshGCM, ui_ex[i]);
			plotVector(meshGCM, ui[i], "GCM_incU_"+i+".dat");
			vi[i] = new SparseVector(dimGCM);
			for(int k=1;k<=dimGCM;k++) {
				vi[i].set(k, Math.log(ui[i].get(k)) / (s[i]*s[i]) );
			}
			plotVector(meshGCM, vi[i], "GCM_incV_"+i+".dat");
		}
		//phi[0] phi[1]
		for(int i=0;i<N-1;i++) {
			phi[i] = new SparseVector(dimGCM);
			for(int k=1;k<=dimGCM;k++) {
				phi[i].set(k, (1.0/h) * ( vi[i].get(k) - vi[i+1].get(k) ));
			}
			plotVector(meshGCM, phi[i], "GCM_phi_"+i+".dat");
		}	
		
		//Tail function  Tail函数（光源最远处）
		//Vector tailT = vi[0];//extractData(meshForward, meshGCM, bkUL);
		//实际上，应该来自从tail重构出来的a(x)
		//Tail改变高度：超过实际高度，GCM方法会降低部分高度，以达到或接近实际高度
		//Tail改变位置（实际(2.0,2.7)）：靠下(2.0,2.6)偏高，靠右(3.0,2.6)偏高，结果显示GCM对a(x)没有任何贡献
		//setMu_a(3.0, 2.6, 0.3, 
		//		2.0, 1); //假设Mu_a高度只有30%；或高度为200%
		setMu_a(1.0, 2.9, 0.3, 
				2.0, 1); //
		setDelta(1.0, 3.5);
		Vector u_tail = solveForwardNeumann(meshForward); //在meshGCM上求解，top边界条件不好确定
		//截取meshForward的部分解到meshGCM上
		u_tail = extractData(meshForward, meshGCM, u_tail);
		Vector tailT = new SparseVector(dimGCM);
		for(int k=1;k<=dimGCM;k++) {
			tailT.set(k, Math.log(u_tail.get(k)) / (s[0]*s[0]) );
		}
		plotVector(meshGCM, tailT, "tailT.dat");

		/*		
		//Check \Laplace{v} + s^2*|\Nabla{v}|^2 = a(x)/s^2
		Vector vecMu_a = discreteFunction(meshGCM,this.mu_a);
		Vector kMu_a = new Vector(dimGCM);
		kMu_a = Vector.axpy(-1.0/(s[0]*s[0]*coefK), vecMu_a, kMu_a);
		plotVector(meshGCM, kMu_a, "kMu_a.dat");

		Vector T_laplace = computeLaplace2D(meshGCM, tailT);
		Vector T_x = computeDerivative(meshGCM, tailT, "x");
		Vector T_y = computeDerivative(meshGCM, tailT, "y");
		Vector v_laplace = T_laplace;
		Vector v_square = Vector.axpy(1.0, 
				Vector.axmy(1.0, T_x, T_x), 
				Vector.axmy(1.0, T_y, T_y)
				);
		plotVector(meshGCM, v_square, "v_square.dat");
		Vector v_equation = Vector.axpy(s[0]*s[0], v_square, v_laplace);
		v_equation = Vector.axpy(-1.0 / (s[0]*s[0]*coefK), vecMu_a, v_equation);
		plotVector(meshGCM, v_equation, "v_equation.dat");
		
		
		//v(x,s[1])
		v_laplace = computeLaplace2D(meshGCM, vi[1]);
		plotVector(meshGCM, v_laplace, "v1_laplace.dat");
		Vector v_x = computeDerivative(meshGCM, vi[1], "x");
		Vector v_y = computeDerivative(meshGCM, vi[1], "y");
		v_square = Vector.axpy(1.0, 
				Vector.axmy(1.0, v_x, v_x), 
				Vector.axmy(1.0, v_y, v_y)
				);
		plotVector(meshGCM, v_square, "v1_square.dat");
		
		Vector phi_laplace = computeLaplace2D(meshGCM, phi[0]);
		Vector phi_x = computeDerivative(meshGCM, phi[0], "x");
		Vector phi_y = computeDerivative(meshGCM, phi[0], "y");
		plotVector(meshGCM, phi_laplace, "phi_laplace.dat");
		
		Vector q_laplace = phi_laplace;
		Vector nablaQ_nablaV = Vector.axpy(1.0,
				Vector.axmy(1.0, phi_x, v_x),
				Vector.axmy(1.0, phi_y, v_y)
				);
		plotVector(meshGCM, nablaQ_nablaV, "nablaQ_nablaV.dat");
		Vector q_equation = new Vector(dimGCM);
		System.out.println(2.0*s[1]*s[1]+"*nablaQ_nablaV");
		q_equation = Vector.axpy(2.0*s[1]*s[1], nablaQ_nablaV, q_laplace);
		System.out.println(4.0*s[1]+"*v_square");
		q_equation = Vector.axpy(4.0*s[1], v_square, q_equation);
		System.out.println(2.0/s[1]+"*v_laplace");
		q_equation = Vector.axpy(2.0/s[1], v_laplace, q_equation);
		//内部为0，边界不为0
		plotVector(meshGCM, q_equation, "q_equation.dat");
		*/
		
		//TODO ???
		setDelta(1.0+(N-1)*h, 3.5);
		solveGCM(meshGCM, N, s, phi, tailT);	
	}
	
	/**
	 * 
	 * @param mesh
	 * @param phi 边界条件
	 * @param tailT
	 */
	/**
	 * 
	 * @param mesh
	 * @param N
	 * @param s
	 * @param phi
	 * @param tailT
	 */
	public  void solveGCM(Mesh mesh, int N, double[]s, Vector[]phi, Vector tailT) {
		int dim = tailT.getDim();
		Vector sumLaplaceQ = new SparseVector(dim);
		Vector sumQ_x= new SparseVector(dim);
		Vector sumQ_y= new SparseVector(dim);
		Vector sumLaplaceQ_real = new SparseVector(dim);
		Vector sumQ_x_real = new SparseVector(dim);
		Vector sumQ_y_real = new SparseVector(dim);
		Vector[] q = new Vector[N];
		Vector T_laplace = Tools.computeLaplace2D(mesh, tailT);
		Vector T_x = Tools.computeDerivative(mesh, tailT, "x");
		Vector T_y = Tools.computeDerivative(mesh, tailT, "y");
		plotVector(mesh, T_laplace, "T_laplace.dat");
		plotVector(mesh, T_x, "T_x.dat");
		plotVector(mesh, T_y, "T_y.dat");
		
		Vector q_0 = new SparseVector(dim);
		for(int i=0;i<N-1;i++) {
			q[i] = solveGCM(mesh,i,
					s[i],s[i+1],
					sumLaplaceQ,sumQ_x,sumQ_y,
					//sumLaplaceQ_real,sumQ_x_real,sumQ_y_real,
					//computeLaplace2D(mesh, sumLaplaceQ_real),
					//computeDerivative(mesh, sumQ_x_real, "x"),
					//computeDerivative(mesh, sumQ_y_real, "y"),
					T_laplace,T_x,T_y,
					new Vector2Function(phi[i]),
					q_0,phi[i]);
		
			q_0 = q[i];
			
			//sum_{j=1}^{N-1} \Laplace{q}
			sumLaplaceQ = SparseVector.axpy(1.0, Tools.computeLaplace2D(mesh, q[i]), sumLaplaceQ);
			sumQ_x = SparseVector.axpy(1.0, Tools.computeDerivative(mesh, q[i], "x"), sumQ_x);
			plotVector(mesh, sumLaplaceQ, "GCM_sumLaplaceQ_"+i+".dat");
			plotVector(mesh, sumQ_x, "GCM_sumQ_x_"+i+".dat");
			plotVector(mesh, sumQ_y, "GCM_sumQ_y_"+i+".dat");
			
			//TODO !!!TEST!!!
//			sumLaplaceQ_real = Vector.axpy(1.0, computeLaplace2D(mesh, phi[i]), sumLaplaceQ_real);
//			sumQ_x_real = Vector.axpy(1.0, computeDerivative(mesh, phi[i], "x"), sumQ_x_real);
//			sumQ_y_real = Vector.axpy(1.0, computeDerivative(mesh, phi[i], "y"), sumQ_y_real);
//			plotVector(mesh, sumLaplaceQ_real, "GCM_sumLaplaceQ_real_"+i+".dat");
//			plotVector(mesh, sumQ_x_real, "GCM_sumQ_x_real_"+i+".dat");
//			plotVector(mesh, sumQ_y_real, "GCM_sumQ_y_real_"+i+".dat");
			
//			sumLaplaceQ_real = Vector.axpy(1.0, phi[i], sumLaplaceQ_real);
//			sumQ_x_real = Vector.axpy(1.0, phi[i], sumQ_x_real);
//			sumQ_y_real = Vector.axpy(1.0, phi[i], sumQ_y_real);
//
//			plotVector(mesh, sumLaplaceQ_real, "GCM_sumLaplaceQ_real_"+i+".dat");
//			plotVector(mesh, sumQ_x_real, "GCM_sumQ_x_real_"+i+".dat");
//			plotVector(mesh, sumQ_y_real, "GCM_sumQ_y_real_"+i+".dat");
			
		}
		
		//用最近处光源重构(注意this.delta的位置)：v_tidle -> u ->　a
		Vector v_tidle = new SparseVector(dim);
		for(int i=0;i<N-1;i++) {
			double h = s[i] - s[i+1];
			v_tidle = SparseVector.axpy(-s[N-1]*s[N-1]*h, q[i], v_tidle);
		}
		plotVector(mesh, v_tidle, "v_tidle_N_befor_add_tail.dat");
		v_tidle = SparseVector.axpy(s[N-1]*s[N-1], tailT, v_tidle);
		plotVector(mesh, v_tidle, "v_tidle_N.dat");
		
		Vector u = new SparseVector(dim);
		for(int i=1;i<=dim;i++) {
			u.set(i, Math.pow(Math.E, v_tidle.get(i)));
		}
		plotVector(mesh, u, "u_N.dat");
		
		Vector a = solveParamInverse(mesh, u);
		plotVector(mesh, a, "alpha_recon_N.dat");
		
		
		//模拟计算出来的phi除了提供边界条件，区域内部的值也有，可以用来重构真实的a(x)
		Vector v_tidle_real = new SparseVector(dim);
		for(int i=0;i<N-1;i++) {
			double h = s[i] - s[i+1];
			v_tidle_real = SparseVector.axpy(-s[N-1]*s[N-1]*h, phi[i], v_tidle_real);
		}
		v_tidle_real = SparseVector.axpy(s[N-1]*s[N-1], tailT, v_tidle_real);
		plotVector(mesh, v_tidle_real, "v_tidle_N_real.dat");
		
		Vector u_real = new SparseVector(dim);
		for(int i=1;i<=dim;i++) {
			u_real.set(i, Math.pow(Math.E, v_tidle_real.get(i)));
		}
		plotVector(mesh, u_real, "u_N_real.dat");
		
		Vector a_real = solveParamInverse(mesh, u_real);
		plotVector(mesh, a_real, "alpha_recon_N_real.dat");
		
	}
	
	public Vector solveGCM(Mesh mesh,int n_1,
			double sn_1, double sn,
			Vector sumLaplaceQ,
			Vector sumQ_x, Vector sumQ_y,
			Vector laplaceT,
			Vector T_x, Vector T_y,
			Function diri,
			Vector q_0,
			Vector phi //for check only
			) {
		
		//Mark border type
		HashMap<NodeType, Function> mapNTF = new HashMap<NodeType, Function>();
		mapNTF.clear();
		mapNTF.put(NodeType.Dirichlet, null);
//		mapNTF.put(NodeType.Dirichlet, new FAbstract("x","y"){
//			@Override
//			public double value(Variable v) {
//				double y = v.get("y");
//				if(Math.abs(y - 3.0) < Constant.eps)
//					return 1.0;
//				else
//					return -1.0;
//			}
//		});
//		
//		//mapNTF.put(NodeType.Robin, null);
		
		double h = sn_1 - sn; //错误的写成：sn-1 - sn，导致结果异常，以后变量名称一定不要这样写了
		Vector A = getGCMCoef(sn_1,sn);
		System.out.println("----------");
		for(int i=1;i<=A.getDim();i++) {
			System.out.print("A"+i+"="+A.get(i)+"  ");
		}
		System.out.println("\n----------");
		int dim = sumLaplaceQ.getDim();
		
		//Right hand side
		Vector f1 = new SparseVector(dim); //zero vector
		//f1 = Vector.axpy(A.get(3)*h, sumLaplaceQ,f1);
		//TODO???
		f1 = SparseVector.axpy(-A.get(3)*h, sumLaplaceQ,f1);
		f1 = SparseVector.axpy(A.get(3), laplaceT,f1);
		
		Vector f2x = new SparseVector(dim); //zero vector
		f2x = SparseVector.axpy(h, sumQ_x,f2x);
		f2x = SparseVector.axpy(-1.0, T_x,f2x);
		f2x = SparseVector.axMuly(A.get(4), f2x, f2x);
		plotVector(mesh, f2x, "GCM_f2x_"+n_1+".dat");
		
		Vector f2y = new SparseVector(dim); //zero vector
		f2y = SparseVector.axpy(h, sumQ_y,f2y);
		f2y = SparseVector.axpy(-1.0, T_y,f2y);
		f2y = SparseVector.axMuly(A.get(4), f2y, f2y);
		plotVector(mesh, f2y, "GCM_f2y_"+n_1+".dat");
		
		Vector f = SparseVector.axpy(1.0, f2x, f2y);
		plotVector(mesh, f, "GCM_f2xpf2y_"+n_1+".dat");
		f = SparseVector.axpy(1.0, f1, f);
		plotVector(mesh, f, "GCM_"+n_1+"_f.dat");
		
		//b1, b2
		Vector q_x_k_1 = null;
		Vector q_y_k_1 = null;
		
		//Nonlinear iteration
		//Vector qn_1 = q_0;
		Vector qn_1 = new SparseVector(dim);
		Vector qn = null;
		for(int iter=1;iter<=6;iter++) {
			
			q_x_k_1 = Tools.computeDerivative(mesh, qn_1, "x");
			q_y_k_1 = Tools.computeDerivative(mesh, qn_1, "y");
			plotVector(mesh, q_x_k_1, "GCM_q_x_km1_"+n_1+"_"+iter+".dat");
			plotVector(mesh, q_y_k_1, "GCM_q_y_km1_"+n_1+"_"+iter+".dat");

			Vector b1 = new SparseVector(dim); //zero vector
			b1 = SparseVector.axpy(A.get(2)*h, sumQ_x,b1);
			b1 = SparseVector.axpy(-A.get(2), T_x,b1);
			b1 = SparseVector.axpy(-A.get(1), q_x_k_1,b1);
			
			Vector b2 = new SparseVector(dim); //zero vector
			b2 = SparseVector.axpy(A.get(2)*h, sumQ_y,b2);
			b2 = SparseVector.axpy(-A.get(2), T_y,b2);
			b2 = SparseVector.axpy(-A.get(1), q_y_k_1,b2);
			
			plotVector(mesh, b1, "GCM_"+n_1+"_b1_"+iter+".dat");
			plotVector(mesh, b2, "GCM_"+n_1+"_b2_"+iter+".dat");
			
			//check lhs with rhs
			if(phi !=  null) {
				Vector phi_laplace = Tools.computeLaplace2D(mesh, phi);
				Vector phi_x = Tools.computeDerivative(mesh, phi, "x");
				Vector phi_y = Tools.computeDerivative(mesh, phi, "y");
				plotVector(mesh, phi_laplace, "GCM_"+n_1+"_phi_laplace_"+iter+".dat");
				plotVector(mesh, phi_x, "GCM_"+n_1+"_phi_x_"+iter+".dat");
				plotVector(mesh, phi_y, "GCM_"+n_1+"_phi_y_"+iter+".dat");
				
				Vector b1r = new SparseVector(dim); //zero vector
				b1r = SparseVector.axpy(A.get(2)*h, sumQ_x, b1r);
				b1r = SparseVector.axpy(-A.get(2),  T_x,    b1r);
				b1r = SparseVector.axpy(-A.get(1),  phi_x,  b1r);
				
				Vector b2r = new SparseVector(dim); //zero vector
				b2r = SparseVector.axpy(A.get(2)*h, sumQ_y, b2r);
				b2r = SparseVector.axpy(-A.get(2),  T_y,    b2r);
				b2r = SparseVector.axpy(-A.get(1),  phi_y,  b2r);

				Vector lhs_x = SparseVector.axMuly(1.0, phi_x, b1r);
				Vector lhs_y = SparseVector.axMuly(1.0, phi_y, b2r);
				Vector lhs = SparseVector.axpy(1.0, lhs_x, lhs_y);
				plotVector(mesh, lhs, "GCM_"+n_1+"_lhsx__lhsy_"+iter+".dat");
				lhs = SparseVector.axpy(1.0, phi_laplace, lhs);
				plotVector(mesh, lhs, "GCM_"+n_1+"_lhs_"+iter+".dat");
				
				//TODO !!!TEST!!!
				//b1 = b1r;
				//b2 = b2r;
			}
			
			WeakFormGCM weakForm = new WeakFormGCM();
			weakForm.setF(new Vector2Function(f));
			weakForm.setParam(
					new FC(-1.0),//注意，有负号!!!
					new FC(0.0), //-\eps{q_n}!!! \eps-0.1
					new Vector2Function(b1),
					new Vector2Function(b2)
				);
			//weakForm.setRobin(new FConstant(0.05),null); 
			
			mesh.clearBorderNodeMark();
			mesh.markBorderNode(mapNTF);

			AssemblerScalar assembler = new AssemblerScalar(mesh, weakForm);
			System.out.println("Begin Assemble...solveGCM");
			Matrix stiff = assembler.getStiffnessMatrix();
			Vector load = assembler.getLoadVector();
			//Dirichlet condition
			assembler.imposeDirichletCondition(diri);
			System.out.println("Assemble done!");

			Solver solver = new Solver();
			qn = solver.solve(stiff, load);
			plotVector(mesh, qn, "GCM_qn_"+n_1+"_"+iter+".dat");
			
			
			Vector q_x_k = Tools.computeDerivative(mesh, qn, "x");
			Vector q_y_k = Tools.computeDerivative(mesh, qn, "y");
			Vector lhs_x = SparseVector.axMuly(1.0, q_x_k, b1);
			Vector lhs_y = SparseVector.axMuly(1.0, q_y_k, b2);
			Vector lhs = SparseVector.axpy(1.0, lhs_x, lhs_y);
			plotVector(mesh, lhs, "GCM_"+n_1+"_lhsx__lhsy_cmp_"+iter+".dat");
			Vector qn_laplace = Tools.computeLaplace2D(mesh, qn);
			lhs = SparseVector.axpy(1.0, qn_laplace, lhs);
			plotVector(mesh, lhs, "GCM_"+n_1+"_lhs_cmp_"+iter+".dat");
			Vector qn_equation = SparseVector.axpy(-1.0, f, lhs);
			plotVector(mesh, qn_equation, "GCM_"+n_1+"_qn_equation_cmp_"+iter+".dat");
			
			
			System.out.println("GCM Iteration "+iter);
			if(iter > 1) {
				double norm = SparseVector.axpy(-1.0, qn_1, qn).norm2();
				if(norm < 0.1) {
					System.out.println("GCM norm < 1.0, STOP");
					break;
				}
			}
			qn_1 = qn;
		}
		return qn;
	}
	
}
