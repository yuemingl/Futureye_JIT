package edu.uta.futureye.application;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;

import edu.uta.futureye.algebra.SolverJBLAS;
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
import edu.uta.futureye.io.MeshReader;
import edu.uta.futureye.io.MeshWriter;
import edu.uta.futureye.lib.assembler.AssemblerScalar;
import edu.uta.futureye.lib.assembler.AssemblerScalarFast;
import edu.uta.futureye.lib.shapefun.SFBilinearLocal2D;
import edu.uta.futureye.lib.shapefun.SFLinearLocal2D;
import edu.uta.futureye.lib.shapefun.SFQuadraticLocal2D;
import edu.uta.futureye.lib.shapefun.SFQuadraticLocal2DFast;
import edu.uta.futureye.lib.weakform.WeakFormL22D;
import edu.uta.futureye.lib.weakform.WeakFormLaplace2D;
import edu.uta.futureye.util.Constant;
import edu.uta.futureye.util.FutureyeException;
import edu.uta.futureye.util.PairDoubleInteger;
import edu.uta.futureye.util.Utils;
import edu.uta.futureye.util.container.ElementList;
import edu.uta.futureye.util.container.NodeList;
import edu.uta.futureye.util.container.ObjList;

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
public class GCMModel {
	public String outputFolder = "";
	
	//Light source
	public Function delta = null;
	public Variable lightSource = null; //light source position
	
	//Inclusion mu_a
	public Function mu_a = null;
	
	//Inclusion 1/(3*mu_s') = 1.0/30.0 ?
	public Function k = new FC(0.02);
	
	public boolean debug = false;
	
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
		delta = new FDelta(this.lightSource,0.01,2e5);
		//测试将dleta函数变得平缓
		//delta = new FDelta(this.lightSource,0.05,2e5);
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
		weakForm.setF(this.k.M(this.delta));
		
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

		SolverJBLAS solver = new SolverJBLAS();
		Vector u = solver.solveDGESV(stiff, load);
		return u;
	}	

	/**
	 * 
	 * @param mesh
	 * @param diri
	 * @param borderType 1: only top boundary is Dirichlet, 2: all boundary are Dirichlet
	 * @return
	 */
	public Vector solveForwardDirichlet(Mesh mesh, Function diri, int borderType) {
		//Mark border type
		HashMap<NodeType, Function> mapNTF = new HashMap<NodeType, Function>();
		mapNTF.clear();
		if(borderType == 1) {
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
		} else {
			mapNTF.put(NodeType.Dirichlet, null);
		}
		mesh.clearBorderNodeMark();
		mesh.markBorderNode(mapNTF);
		
		WeakFormLaplace2D weakForm = new WeakFormLaplace2D();
		
		//Right hand side
		weakForm.setF(this.k.M(this.delta));
		//rhs=0 与非零没有差别
		//weakForm.setF(FC.c0);

		weakForm.setParam(
				this.k, this.mu_a, null, this.k //d==k,q=0 (即：u_n + u =0)
			);
		
		//bugfix 2011-5-7两种方式结果不一样？
		//Assembler assembler = new AssemblerScalarFast(mesh, weakForm);
		Assembler assembler = new AssemblerScalar(mesh, weakForm);
		System.out.println("Begin Assemble...solveForwardDirichlet");
		assembler.assemble();
		Matrix stiff = assembler.getStiffnessMatrix();
		Vector load = assembler.getLoadVector();
		//Dirichlet condition
		assembler.imposeDirichletCondition(diri);
		
		System.out.println("Assemble done!");

		SolverJBLAS solver = new SolverJBLAS();
		Vector u = solver.solveDGESV(stiff, load);
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
		weakFormL2.setF(this.k.M(this.delta));
		//weakFormL2.setF(FC.c0); //bugfix 2011-5-7 不能用0
		
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
		
		SolverJBLAS solver = new SolverJBLAS();
		Vector u = solver.solveDGESV(stiff, load);
		return u;
	}
	
	public Vector solveParamInverse2(Mesh mesh, Vector U, Function f) {
		HashMap<NodeType, Function> mapNTF2 = new HashMap<NodeType, Function>();
		mapNTF2.put(NodeType.Dirichlet, null);
		mesh.clearBorderNodeMark();
		mesh.markBorderNode(mapNTF2);
		
		//Weak form
		WeakFormL22D weakFormL2 = new WeakFormL22D();
		//Right hand side
		weakFormL2.setF(f);
		
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
		
		SolverJBLAS solver = new SolverJBLAS();
		Vector u = solver.solveDGESV(stiff, load);
		return u;
	}
	/**
	 * 
	 * @param mesh
	 * @param a_m
	 * @param a_m_1
	 * @param u_m_1
	 * @param iterNum
	 * @param boundary 边界条件：边界条件为0的假设在back reflect时有问题，不为0
	 * @return
	 */
	public Vector solveEnhance(Mesh mesh, 
			Vector a_m, Vector a_m_1, Vector u_m_1,
			int iterNum, Vector boundary
			) {
		//Mark border type
		HashMap<NodeType, Function> mapNTF = new HashMap<NodeType, Function>();
		mapNTF.clear();
//		mapNTF.put(NodeType.Robin, null);
		mapNTF.put(NodeType.Dirichlet, null);
		
		mesh.clearBorderNodeMark();
		mesh.markBorderNode(mapNTF);
		
		WeakFormLaplace2D weakForm = new WeakFormLaplace2D();
		
		//Right hand side
		Function rhs = null;
		Function diff_am = new Vector2Function(FMath.axpy(-1.0, a_m_1, a_m));
		this.plotFunction(mesh, diff_am, "enhance_diff_am"+iterNum+".dat");
		if(iterNum>=2) {
			Function lamd_exp = new FC(Math.PI*Math.PI*Math.pow(Math.E, iterNum-1)).
								M(diff_am);
			//lamd_exp = new FConstant(0.8);
			this.plotFunction(mesh, lamd_exp, "enhance_rhs_lamd_exp"+iterNum+".dat");
			Function lamd = FMath.pow(new FC(Math.E), lamd_exp).M(new FC(Math.pow(1.05, -iterNum)));
			this.plotFunction(mesh, lamd, "enhance_rhs_lamd"+iterNum+".dat");
			rhs = lamd.M(diff_am).M(new Vector2Function(u_m_1));
			this.plotFunction(mesh, rhs, "enhance_rhs"+iterNum+".dat");
		} else {
			rhs = diff_am.M(new Vector2Function(u_m_1)).M(1.5);//enhance!!!
			this.plotFunction(mesh, rhs, "enhance_rhs"+iterNum+".dat");
		}
		weakForm.setF(FC.c0.S(rhs)); //-1
		
//		weakForm.setParam(
//			new FC(-1.0), //2011-5-7 负系数
//			new FC(0.0).S(new Vector2Function(a_m)),
//			null,null
//		);
		//2011-5-7 OK
		weakForm.setParam(
				this.k, new Vector2Function(a_m), null, this.k
		);
		
		//bugfix 2011-5-7两种方式结果不一样？
		//Assembler assembler = new AssemblerScalarFast(mesh, weakForm);
		Assembler assembler = new AssemblerScalar(mesh, weakForm);
		System.out.println("Begin Assemble...");
		assembler.assemble();
		Matrix stiff = assembler.getStiffnessMatrix();
		Vector load = assembler.getLoadVector();
		if(boundary != null)
			assembler.imposeDirichletCondition(new Vector2Function(boundary));
		else
			assembler.imposeDirichletCondition(new FC(0.0));
		System.out.println("Assemble done!");

		SolverJBLAS solver = new SolverJBLAS();
		Vector u = solver.solveDGESV(stiff, load);
		
	    Vector alpha_wL = solveParamInverse2(mesh,u,FC.c0.S(rhs));
	    plotVector(mesh, alpha_wL, "enhance_alpha_wL"+iterNum+".dat");

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
	 * 会改变网格的边界条件！！！
	 * @param mesh
	 * @param bkUL
	 * @param incUL
	 * @return
	 */
	public Vector computeTailLeftLightSource(Mesh mesh,Vector bkUL, Vector incUL) {
		//Mark border nodes
		HashMap<NodeType, Function> mapNTF = new HashMap<NodeType, Function>();
		mapNTF.clear();
		mapNTF.put(NodeType.Robin, null);
		mesh.clearBorderNodeMark();
		mesh.markBorderNode(mapNTF);
		
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

	
	/**
	 * 依赖于mesh的大小！！！
	 * 会改变网格的边界条件！！！
	 * @param mesh
	 * @param bkUL
	 * @param incUL
	 * @return
	 */
	public Vector computeTailLeftLightSource2(Mesh mesh,Vector bkUL, Vector incUL) {
		//Mark border nodes
		HashMap<NodeType, Function> mapNTF = new HashMap<NodeType, Function>();
		mapNTF.clear();
		mapNTF.put(NodeType.Robin, null);
		mesh.clearBorderNodeMark();
		mesh.markBorderNode(mapNTF);
		
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
		//Mark border nodes
		HashMap<NodeType, Function> mapNTF = new HashMap<NodeType, Function>();
		mapNTF.clear();
		mapNTF.put(NodeType.Robin, null);
		mesh.clearBorderNodeMark();
		mesh.markBorderNode(mapNTF);

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
	
	/**
	 * 为网格单元赋值自由度
	 * @param elementType
	 * @param mesh
	 */
	public void assignDOF(int elementType, Mesh mesh) {
		ScalarShapeFunction[] shapeFun = null;
		ScalarShapeFunction[] shapeFunRect = null;
		if(elementType == 1) {
			mesh.computeNodeBelongsToElements();
			mesh.computeNeighborNodes();
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
			}
			mesh.computeNodeBelongsToElements();
			mesh.computeNeighborNodes();
			//Assign degree of freedom to element
			shapeFun = new SFQuadraticLocal2DFast[6];
			for(int i=1;i<=6;i++)
				shapeFun[i-1] = new SFQuadraticLocal2DFast(i);
//			shapeFun = new SFQuadraticLocal2D[6];
//			for(int i=1;i<=6;i++)
//				shapeFun[i-1] = new SFQuadraticLocal2D(i);

			//Assign shape function to DOF
			for(int i=1;i<=mesh.getElementList().size();i++) {
				Element e = mesh.getElementList().at(i);
				for(int j=1;j<=e.nodes.size();j++) {
					DOF dof = new DOF(j,e.nodes.at(j).globalIndex,shapeFun[j-1]);
					e.addNodeDOF(j, dof);
				}
			}
		} else {
			throw new FutureyeException("Error: elementType parameter!");
		}		
	}
	
	/**
	 * 入口函数
	 * 
	 * @param gridFileForward
	 * @param elementType
	 * @param gridFileGCM
	 * @param outputFolder
	 */
	public void run(int elementType, String gridFileForward, String gridFileGCM, 
			String outputFolder) {
		MeshReader reader = null;
		reader = new MeshReader(gridFileForward);
		Mesh meshForward = reader.read2DMesh();
		reader = new MeshReader(gridFileGCM);
		Mesh meshGCM = reader.read2DMesh();
		this.outputFolder = outputFolder;
		this.assignDOF(elementType, meshForward);
		this.assignDOF(elementType, meshGCM);

		NodeList list = meshGCM.getNodeList();
		int nNode = list.size();
		double deltaY = 3.5; //光源的y坐标
		double mu_aY = 2.7; //inclustion的y坐标
		double mu_aMax = 2.0;
		
		//----------------------Begin Left light source ----------------------------
		setDelta(1.0, deltaY);
		
		//Solve background forward problem
		setMu_a(0.0, 0.0, 0.0, 
				0.1, 1);
		Vector bkUL = solveForwardNeumann(meshForward);
		plotVector(meshForward, bkUL, "bkUL.dat");
		
		//Solve forward problem with inclusion
		setMu_a(2.0, mu_aY, 0.3, 
				mu_aMax, 1);
		plotFunction(meshForward, this.mu_a, "alpha_real.dat");
		Vector incUL = solveForwardNeumann(meshForward);
		plotVector(meshForward, incUL, "incUL.dat");
		incUL = this.extractData(meshForward, meshGCM, incUL);
		bkUL = this.extractData(meshForward, meshGCM, bkUL);
		//会改变网格的边界条件!!!
	    Vector tailUL = computeTailLeftLightSource(meshGCM, bkUL, incUL);
	    plotVector(meshGCM, tailUL, "tailUL.dat");
	    
		//Recovery parameter mu_a from solution
	    Vector alpha_real_cmpt = solveParamInverse(meshGCM,incUL);
	    plotVector(meshGCM, alpha_real_cmpt, "alpha_real_cmpt.dat");

	    //Recovery parameter mu_a from tail
	    Vector tailUL_noise = addNoise(tailUL,0.0);
	    Vector alphaL = solveParamInverse(meshGCM,tailUL_noise);
	    //Cut noise
	    for(int i=1;i<=nNode;i++) {
	    	Node node = list.at(i);
	    	if(node.coord(1) <1.3 || node.coord(1)>3.7 || node.coord(2)<2.0) {
		    //if(node.coord(1) <1.9 || node.coord(1)>2.8 || node.coord(2)<1.5) {
	    		alphaL.set(node.globalIndex, 0.1);
	    	}
		}
	    plotVector(meshGCM, alphaL, "alphaL.dat");

	    
		//----------------------Begin Right light source ---------------------
		setDelta(4.0, deltaY);
		
		setMu_a(0.0, 0.0, 0.0, 
				0.1, 1);
		Vector bkUR = solveForwardNeumann(meshForward);
		plotVector(meshForward, bkUR, "bkUR.dat");
		
		setMu_a(2.0, mu_aY, 0.3, 
				mu_aMax, 1);
		Vector incUR = solveForwardNeumann(meshForward);
	    plotVector(meshForward, incUR, "incUR.dat");
		incUR = this.extractData(meshForward, meshGCM, incUR);
		bkUR = this.extractData(meshForward, meshGCM, bkUR);
		//会改变网格的边界条件!!!
	    Vector tailUR = computeTailRightLightSource(meshGCM, bkUR, incUR);
	    plotVector(meshGCM, tailUR, "tailUR.dat");
	    
	    //Recovery parameter mu_a from tail
	    Vector alphaR = solveParamInverse(meshGCM,tailUR);
	    //Cut noise
	    for(int i=1;i<=nNode;i++) {
	    	Node node = list.at(i);
	    	if(node.coord(1) <1.3 || node.coord(1)>3.7 || node.coord(2)<2.0) {
		    //if(node.coord(1) <1.9 || node.coord(1)>2.8 || node.coord(2)<1.5) {
	    		alphaR.set(node.globalIndex, 0.1);
	    	}
		}
	    plotVector(meshGCM, alphaR, "alphaR.dat");
	    
	    //Smooth alphaL and alphaR...
	    alphaL = Utils.gaussSmooth(meshGCM, alphaL, 1, 0.5);
	    alphaL = Utils.gaussSmooth(meshGCM, alphaL, 1, 0.5);
	    alphaL = Utils.gaussSmooth(meshGCM, alphaL, 1, 0.5);
	    alphaR = Utils.gaussSmooth(meshGCM, alphaR, 1, 0.5);
	    alphaR = Utils.gaussSmooth(meshGCM, alphaR, 1, 0.5);
	    alphaR = Utils.gaussSmooth(meshGCM, alphaR, 1, 0.5);

	    Vector alpha_avg = new SparseVector(nNode);
	    for(int i=1;i<=nNode;i++) {
	    	alpha_avg.set(i, (alphaL.get(i)+alphaR.get(i))/2.0);
	    }
	    plotVector(meshGCM, alpha_avg, "alpha_avg.dat");
	    
	    //Smooth alpha_avg...
	    alpha_avg = Utils.gaussSmooth(meshGCM, alpha_avg, 2, 0.5);
	    alpha_avg = Utils.gaussSmooth(meshGCM, alpha_avg, 2, 0.5);
	    alpha_avg = Utils.gaussSmooth(meshGCM, alpha_avg, 2, 0.5);
	    Vector alpha_avg_smooth = Utils.gaussSmooth(meshGCM, alpha_avg, 2, 0.5);
	    plotVector(meshGCM, alpha_avg_smooth, "alpha_avg_smooth.dat");
	    //Cut 80%
	    Double max = FMath.max(alpha_avg_smooth);
	    for(int i=1;i<=nNode;i++) {
	    	if(alpha_avg_smooth.get(i) < 0.8*max || alpha_avg_smooth.get(i)<0.1 ||
	    			!list.at(i).isInnerNode())
	    		alpha_avg_smooth.set(i, 0.1);
	    	else
	    		alpha_avg_smooth.add(i, 0.0);
	    }
	    plotVector(meshGCM, alpha_avg_smooth, "alpha_avg_smooth_cut.dat");
	    
	    
	    //-------------------- Begin Nonlinear Enhancement -----------------------

	    setDelta(1.0, deltaY);
	    //方法一：
	    /*
	    Vector alpha_m1 = alpha_avg_smooth;
		//边界条件需要满足测量值，而不是Neumann边界，但是结果很差（bugfix：结果很好 solveForwardDirichlet()：AssemblerScalarFast->AssemblerScalar）
		mu_a = new Vector2Function(alpha_m1);
		//setMu_a(2.0, mu_aY, 0.3, 
		//		mu_aMax, 1);
		//alpha_m1 = alpha_real_cmpt;
		Function diri = new Vector2Function(incUL);
		Vector um1 = solveForwardDirichlet(meshGCM,diri,2);
		
		//在meshForward上使用Neumann边界，然后截取的meshGCM上反而是正确的（bugfix：结果待定）
		//mu_a = new Vector2Function(this.extendData(meshGCM, meshForward, alpha_m1, 0.1));
		//this.plotFunction(meshForward, mu_a, "mu_a_extend.dat");
		//Vector um1= this.solveForwardNeumann(meshForward);
		//um1 = this.extractData(meshForward, meshGCM, um1);
	    plotVector(meshGCM, um1, "enhance_um1.dat");
	    
	    Vector incUL_subtract_bkUL = FMath.axpy(-1.0, bkUL, incUL);
	    plotVector(meshGCM, incUL_subtract_bkUL, "incUL_subtract_bkUL.dat");

		for(int nit=2;nit<=3;nit++) {
		    Vector alpha_m2 = solveParamInverse(meshGCM,um1);
			//可以考虑cut，问题：alpha_m2==alpha_m1，会导致无法enhance
			plotVector(meshGCM, alpha_m2, "enhance_alpha_m"+nit+".dat");

		    Vector wL1 = solveEnhance(meshGCM,alpha_m2,alpha_m1,um1,nit,
		    		incUL_subtract_bkUL);
		    //Vector wL1 = solveEnhance(meshGCM,alpha_m2,alpha_m1,um1,nit,
		    //		null);
		    plotVector(meshGCM, wL1, "enhance_wL"+nit+".dat");
		    
		    Vector um2 = FMath.axpy(1.0, um1, wL1);
		    plotVector(meshGCM, um2, "enhance_um"+nit+".dat");
		    
		    um1 = um2;
		    alpha_m1 = alpha_m2;
		}
		*/
	    
	    //方法二：从背景开始
	    Vector um1 = bkUL;
		//初值使用背景
		Vector alpha_m1 = new SparseVector(nNode,0.1);
	    plotVector(meshGCM, alpha_m1, "enhance_alpha_m01.dat");
		Vector alpha_m2 = alpha_avg_smooth;
	    plotVector(meshGCM, alpha_m2, "enhance_alpha_m02.dat");
	    Vector incUL_subtract_bkUL = FMath.axpy(-1.0, bkUL, incUL);
	    plotVector(meshGCM, incUL_subtract_bkUL, "incUL_subtract_bkUL.dat");
	    
		mu_a = new Vector2Function(alpha_m1);
		Function diri = new Vector2Function(incUL);
		Vector um01 = solveForwardDirichlet(meshGCM,diri,2);
		Vector um2 = null;
		for(int nit=1;nit<=3;nit++) {
		    //Vector wL1 = solveEnhance(meshGCM,alpha_m2,alpha_m1,um1,nit,
		    //		incUL_subtract_bkUL);
		    Vector wL1 = solveEnhance(meshGCM,alpha_m2,alpha_m1,um1,nit,
		    		null);
		    plotVector(meshGCM, wL1, "enhance_wL"+nit+".dat");
		    
		    if(nit==1)
		    	um2 = FMath.axpy(1.0, um01, wL1);
		    else
		    	um2 = FMath.axpy(1.0, um1, wL1);

		    plotVector(meshGCM, um2, "enhance_um"+nit+".dat");
		    
		    alpha_m2 = solveParamInverse(meshGCM,um2);
			//可以考虑cut，问题：alpha_m2==alpha_m1，会导致无法enhance
			plotVector(meshGCM, alpha_m2, "enhance_alpha_m"+nit+".dat");
		    
		    um1 = um2;
		    alpha_m1 = alpha_m2;
		}
		
		//-----------------Begin Global Convergence Method-------------------------
		int N=15;
		double[]s = new double[N];
		double h = 0.06;
		s[0] = 3.0;
		//由远到近，s递减
		for(int i=1; i<N; i++)
			s[i] = s[0] - i*h;
		
		Vector tailT = new SparseVector(um1.getDim());
		//tail函数（最远处）
		for(int i=1;i<=tailT.getDim();i++) {
			tailT.set(i, Math.log(um1.get(i))/(s[0]*s[0]));
		}
		this.plotVector(meshGCM, tailT, "tailT.dat");
		
		//Simulate border measurement data
		Vector[] ui = new Vector[N];
		Vector[] phi = new Vector[N];
		setMu_a(2.0, mu_aY, 0.3, 
				mu_aMax, 1);
		//mu_a = new VectorBasedFunction(alpha_m1);
		for(int i=0;i<N;i++) {
			setDelta(1.0+i*h, deltaY);
			ui[i] = solveForwardNeumann(meshForward);
			ui[i] = this.extractData(meshForward, meshGCM, ui[i]);
			plotVector(meshGCM, ui[i], "GCM_incU_"+i+".dat");
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
			plotVector(meshGCM, phi[i], "GCM_phi_"+i+".dat");
		}
		
		//模拟计算出来的phi除了提供边界条件，区域内部的值也有，可以用来重构真实的a(x)
		setDelta(1.0+(N-1)*h, deltaY);
		solveGCM(meshGCM, N, s, phi, tailT);
		
	}
	
	
	/////////////////////////////Adaptive Mesh Refinement/////////////////////////////
	
	
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
			mesh.computeNodeBelongsToElements();
			mesh.computeNeighborNodes();
			mesh.computeNeighborElements();
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
		plotVector(mesh, bkUL, "bkUL.dat");
		
		//Solve forward problem with inclusion
		setMu_a(2.5, 2.6, 0.3, 
				0.2, 1);
		plotFunction(mesh, this.mu_a, "alpha_real.dat");
		Vector incUL = solveForwardNeumann(mesh);
		plotVector(mesh, incUL, "incUL.dat");
		//中间不要加入其他代码！
	    Vector tailUL = computeTailLeftLightSource(mesh, bkUL, incUL);
	    plotVector(mesh, tailUL, "tailUL.dat");
	    plotTopBorder(mesh,incUL, "incUL_topLine.dat");
	    
		//Recovery parameter mu_a from solution
	    Vector alpha_real_cmpt = solveParamInverse(mesh,incUL);
	    plotVector(mesh, alpha_real_cmpt, "alpha_real_cmpt.dat");

	    //Recovery parameter mu_a from tail
	    Vector tailUL_noise = addNoise(tailUL,0.0);
	    Vector alphaL = solveParamInverse(mesh,tailUL_noise);
	    plotVector(mesh, alphaL, "alphaL.dat");
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
	    plotVector(mesh, alphaL, "alphaL_cut.dat");
	    alphaL = Utils.gaussSmooth(mesh, alphaL, 1, 0.8);
	    alphaL = Utils.gaussSmooth(mesh, alphaL, 1, 0.8);
	    alphaL = Utils.gaussSmooth(mesh, alphaL, 1, 0.8);
	    alphaL = Utils.gaussSmooth(mesh, alphaL, 1, 0.8);
	    plotVector(mesh, alphaL, "alphaL_cut_smooth.dat");

	    	    
		//----------------------Begin Right light source ---------------------
		setDelta(4.0, 2.8);
		
		setMu_a(0.0, 0.0, 0.0, 
				0.1, 1);
		Vector bkUR = solveForwardNeumann(mesh);
		plotVector(mesh, bkUR, "bkUR.dat");
		
		setMu_a(2.5, 2.75, 0.15, 
				2.0, 1);
		Vector incUR = solveForwardNeumann(mesh);
	    plotVector(mesh, incUR, "incUR.dat");
	    Vector tailUR = computeTailRightLightSource(mesh, bkUR, incUR);
	    plotVector(mesh, tailUR, "tailUR.dat");
	    
	    Vector alphaR = solveParamInverse(mesh,tailUR);
	    plotVector(mesh, alphaR, "alphaR.dat");
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
	    plotVector(mesh, alphaR, "alphaR_cut.dat");
	    alphaR = Utils.gaussSmooth(mesh, alphaR, 1, 0.8);
	    alphaR = Utils.gaussSmooth(mesh, alphaR, 1, 0.8);
	    alphaR = Utils.gaussSmooth(mesh, alphaR, 1, 0.8);
	    alphaR = Utils.gaussSmooth(mesh, alphaR, 1, 0.8);
	    plotVector(mesh, alphaR, "alphaR_cut_smooth.dat");
	    
	    Vector alpha_avg = new SparseVector(nNode);
	    for(int i=1;i<=nNode;i++) {
	    	alpha_avg.set(i, (alphaL.get(i)+alphaR.get(i))/2.0);
	    }
	    plotVector(mesh, alpha_avg, "alpha_avg.dat");
	    
	    return alpha_avg;
	}
	
	
	/**
	 * 自适应入口函数
	 * 
	 * @param elementType
	 * @param gridID
	 * @param outputFolder
	 */
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
	    plotVector(mesh, alpha_avg_smooth, "alpha_avg_smooth.dat");
		Vector alpha_avg_smooth_no_noise = alpha_avg_smooth.copy();
	    Double max = alpha_avg_smooth_no_noise.normInf();
	    for(int i=1;i<=alpha_avg_smooth_no_noise.getDim();i++) {
	    	if(alpha_avg_smooth_no_noise.get(i) < 0.7*max)
	    		alpha_avg_smooth_no_noise.set(i, 0.1);
	    }
	    plotVector(mesh, alpha_avg_smooth_no_noise, "alpha_avg_smooth_no_noise.dat");		
		
		//Adaptive refinement 1 -> 输出文件夹：outputFolder_adaptive
		ElementList eToRefine = computeRefineElement(mesh, alpha_avg_smooth, 0.06);
		System.out.println("Before refine: Element="+mesh.getElementList().size()+", Node="+mesh.getNodeList().size());
		Refiner.refineOnce(mesh, eToRefine);
		System.out.println("After refine: Element="+mesh.getElementList().size()+", Node="+mesh.getNodeList().size());
		plotFunction(mesh, this.mu_a, "alpha_real_refine.dat");
	
		alpha_avg = solveAdaptive(mesh, elementType, outputFolder+"\\adaptive");
		alpha_avg_smooth = Utils.gaussSmooth(mesh, alpha_avg, 1, 0.5);
		alpha_avg_smooth = Utils.gaussSmooth(mesh, alpha_avg, 1, 0.5);
		alpha_avg_smooth = Utils.gaussSmooth(mesh, alpha_avg, 1, 0.5);
		alpha_avg_smooth = Utils.gaussSmooth(mesh, alpha_avg, 1, 0.5);
	    alpha_avg_smooth = Utils.gaussSmooth(mesh, alpha_avg, 1, 0.5);
	    alpha_avg_smooth = Utils.gaussSmooth(mesh, alpha_avg, 1, 0.5);
	    plotVector(mesh, alpha_avg_smooth, "alpha_avg_smooth.dat");
	    alpha_avg_smooth_no_noise = alpha_avg_smooth.copy();
	    max = alpha_avg_smooth_no_noise.normInf();
	    for(int i=1;i<=alpha_avg_smooth_no_noise.getDim();i++) {
	    	if(alpha_avg_smooth_no_noise.get(i) < 0.7*max)
	    		alpha_avg_smooth_no_noise.set(i, 0.1);
	    }
	    plotVector(mesh, alpha_avg_smooth_no_noise, "alpha_avg_smooth_no_noise.dat");		
		
		//Adaptive refinement 2 -> 输出文件夹：outputFolder_adaptive2
		eToRefine.clear();
		eToRefine = computeRefineElement(mesh, alpha_avg_smooth, 0.10);
		System.out.println("Before refine: Element="+mesh.getElementList().size()+", Node="+mesh.getNodeList().size());
		Refiner.refineOnce(mesh, eToRefine);
		System.out.println("After refine: Element="+mesh.getElementList().size()+", Node="+	mesh.getNodeList().size());
		plotFunction(mesh, this.mu_a, "alpha_real_refine.dat");
		
		alpha_avg = solveAdaptive(mesh, elementType, outputFolder+"\\adaptive2");
	    alpha_avg_smooth = Utils.gaussSmooth(mesh, alpha_avg, 1, 0.5);
	    alpha_avg_smooth = Utils.gaussSmooth(mesh, alpha_avg, 1, 0.5);
	    alpha_avg_smooth = Utils.gaussSmooth(mesh, alpha_avg, 1, 0.5);
	    plotVector(mesh, alpha_avg_smooth, "alpha_avg_smooth.dat");
	    alpha_avg_smooth_no_noise = alpha_avg_smooth.copy();
	    max = alpha_avg_smooth_no_noise.normInf();
	    for(int i=1;i<=alpha_avg_smooth_no_noise.getDim();i++) {
	    	if(alpha_avg_smooth_no_noise.get(i) < 0.7*max)
	    		alpha_avg_smooth_no_noise.set(i, 0.1);
	    }
	    plotVector(mesh, alpha_avg_smooth_no_noise, "alpha_avg_smooth_no_noise.dat");		
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
	
	///////////////////////////////////GCM////////////////////////////////////////////////
	
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
	
	public Vector extendData(Mesh meshFrom, Mesh meshTo, Vector u, double deaultValue) {
		NodeList nodeTo = meshTo.getNodeList();
		int dimTo = nodeTo.size();
		Vector rlt = new SparseVector(dimTo);
		for(int i=1;i<=nodeTo.size();i++) {
			Node node = meshFrom.containNode(nodeTo.at(i));
			if(node != null) {
				rlt.set(i, u.get(node.globalIndex));
			} else {
				rlt.set(i, deaultValue);
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
	 * 观察Tail的影响对重构结果有多大？是否像之前的测试结果，影响很大？ (2011-5-10 yes)
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
		this.assignDOF(elementType, meshForward);
		this.assignDOF(elementType, meshGCM);

		double coefK = 0.02;
		double mu_aX = 3.0; //inclustion的x坐标
		double mu_aY = 2.6; //inclustion的y坐标
		double deltaX = 0.5; //光源的x坐标
		

		
		setDelta(deltaX, 3.5);
		//setDelta(0.0, 4.0);
		plotFunction(meshForward, this.delta, "delta.dat");
		
		//Solve background forward problem
		setMu_a(0.0, 0.0, 0.0, 
				0.1, 1);
		Vector bkUL = solveForwardNeumann(meshForward);
		plotVector(meshForward, bkUL, "bkUL.dat");
		
		//Solve forward problem with inclusion
		//实际的a(x)值
		setMu_a(mu_aX, mu_aY, 0.3, 
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
		
		//-k*\Delta{u} + mu_a(x)*u  ( = 0 )检验是否为零
		Vector vMu_a = discreteFunction(meshForward,this.mu_a);
		Vector eq = FMath.axMuly(1.0, vMu_a, incUL);
		plotVector(meshForward, eq, "incUL_Mu_amU.dat");//Mu_a*U
		eq = FMath.axpy(-coefK, incUL_laplace, eq);
		plotVector(meshForward, eq, "incUL_equation.dat");//-k\Delta{u} + mu_a(x)*u


		//---------------Global Convergence Method---------------
		//注意，需要修改delta
		//delta = new FDelta(this.lightSource,0.01,2e5);

		int N=30; //Total number of positions of moving light source x_0
		double[]s = new double[N];
		double h = 0.1; //Step size h, //2011-5-6: h=0.4的时候phi所包含的Inclusion比较强，重构出来的a(x)要高于h=0.1的时候
		s[0] = 4.0; //Farthest position of light source: s[0]>=1?  s[0]=1 NO; s[0]=2,4 OK
		for(int i=1; i<N; i++)
			s[i] = s[0] - i*h; //由远到近，s递减
		
		//Simulated boundary measurement data
		//模拟计算边界测量数据
		int dimGCM = meshGCM.getNodeList().size();
		Vector[] ui_ex = new Vector[N];
		Vector[] ui = new Vector[N];
		Vector[] vi = new Vector[N];
		Vector[] phi = new Vector[N];
		
		//ui[0] ui[1] ui[2] ... ui[N-1]
		//vi[0] vi[1] vi[2] ... vi[N-1], *** vi := ln(ui)/(s^2) ***
		for(int i=0;i<N;i++) {
			//实际的a(x)值（用来模拟计算边界phi）
			setMu_a(mu_aX, mu_aY, 0.3, 
					2.0, 1);
			setDelta(deltaX+i*h, 3.5);
			
			ui_ex[i] = solveForwardNeumann(meshForward); //在meshGCM上求解，top边界条件不好确定
			plotVector(meshForward, ui_ex[i], String.format("GCM_incU_ex_%02d.dat",i));
			
			//截取meshForward的部分解到meshGCM上
			ui[i] = extractData(meshForward, meshGCM, ui_ex[i]);
			plotVector(meshGCM, ui[i], String.format("GCM_incU_%02d.dat",i));
			vi[i] = new SparseVector(dimGCM);
			for(int k=1;k<=dimGCM;k++) {
				vi[i].set(k, Math.log(ui[i].get(k)) / (s[i]*s[i]) );
			}
			plotVector(meshGCM, vi[i], String.format("GCM_incV_%02d.dat",i));
		}
		//Boundary function on \partial\Omega (实际上模拟的phi[i]在整个区域上都有值)
		//phi[0] phi[1] ... phi[N-2], *** phi := (vi[i]-vi[i+1])/h ***
		for(int i=0;i<N-1;i++) {
			phi[i] = new SparseVector(dimGCM);
			for(int k=1;k<=dimGCM;k++) {
				phi[i].set(k, (1.0/h) * ( vi[i].get(k) - vi[i+1].get(k) ));
			}
			plotVector(meshGCM, phi[i], String.format("GCM_phi_%02d.dat",i));
		}
		
		//Compute Tail function  计算Tail函数（光源最远处）
		
		//注意：
		//Tail改变高度：超过实际高度，GCM方法会降低部分高度，以达到或接近实际高度
		//Tail改变位置（实际(2.0,2.7)）：靠下(2.0,2.6)偏高，靠右(3.0,2.6)偏高，结果显示GCM对a(x)没有任何贡献
		//setMu_a(3.0, mu_aY, 0.3, 
		//		2.0, 1); //假设Mu_a高度只有30%；或高度为200%
		
		//Test 1
		//从Tail重构出来的a(x)，这里为了计算Tail，设置一个不准确的a(x)
		//setMu_a(1.0, mu_aY, 0.3, 
		//		1.0, 1); //Tail 的Inclusion x坐标在1.0，实际Inclusion x坐标在2.0，但边界条件是在2.0，测试CGM的重构结果
		setMu_a(mu_aX, mu_aY, 0.3, 
				2.0, 1); //真实Tail，q_n应该与真解相同 (经验证，非常一致2011-5-10)
		
		setDelta(deltaX, 3.5);
		Vector u_tail = solveForwardNeumann(meshForward); //在meshGCM上求解，top边界条件不好确定
		//截取meshForward的部分解到meshGCM上
		u_tail = extractData(meshForward, meshGCM, u_tail);
		Vector tailT = new SparseVector(dimGCM);
		for(int k=1;k<=dimGCM;k++) {
			tailT.set(k, Math.log(u_tail.get(k)) / (s[0]*s[0]) );
		}
		plotVector(meshGCM, tailT, "tailT.dat");

		/**
		//Check \Laplace{v} + s^2*|\Nabla{v}|^2 = a(x)/s^2
		Vector vecMu_a = discreteFunction(meshGCM,this.mu_a);
		Vector kMu_a = new Vector(dimGCM);
		kMu_a = Vector.axpy(-1.0/(s[0]*s[0]*coefK), vecMu_a, kMu_a);
		plotVector(meshGCM, kMu_a, "kMu_a.dat");

		Vector T_laplace = computeLaplace2D(meshGCM, tailT);
		Vector T_x = computeDerivative(meshGCM, tailT, "x");
		Vector T_y = computeDerivative(meshGCM, tailT, "y");
		//Check equation of v(x,s[0]) (v == tailT)
		Vector v_laplace = T_laplace;
		Vector v_square = Vector.axpy(1.0, 
				Vector.axmy(1.0, T_x, T_x), 
				Vector.axmy(1.0, T_y, T_y)
				);
		plotVector(meshGCM, v_square, "v_square.dat");
		Vector v_equation = Vector.axpy(s[0]*s[0], v_square, v_laplace);
		v_equation = Vector.axpy(-1.0 / (s[0]*s[0]*coefK), vecMu_a, v_equation);
		plotVector(meshGCM, v_equation, "v_equation.dat");
		
		
		//Check equation of v(x,s[1])
		v_laplace = computeLaplace2D(meshGCM, vi[1]);
		plotVector(meshGCM, v_laplace, "v1_laplace.dat");
		Vector v_x = computeDerivative(meshGCM, vi[1], "x");
		Vector v_y = computeDerivative(meshGCM, vi[1], "y");
		v_square = Vector.axpy(1.0, 
				Vector.axmy(1.0, v_x, v_x), 
				Vector.axmy(1.0, v_y, v_y)
				);
		plotVector(meshGCM, v_square, "v1_square.dat");
		
		//Check equation of q_1
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
		
		//注意光源的位置，定位在最近位置（s最小值）
		setDelta(deltaX+(N-1)*h, 3.5);
		solveGCM(meshGCM, N, s, phi, tailT);	
	}
	

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
//		Vector sumLaplaceQ_real = new SparseVector(dim);
//		Vector sumQ_x_real = new SparseVector(dim);
//		Vector sumQ_y_real = new SparseVector(dim);
		Vector[] q = new Vector[N];
		Vector T_laplace = Tools.computeLaplace2D(mesh, tailT);
		Vector T_x = Tools.computeDerivative(mesh, tailT, "x");
		Vector T_y = Tools.computeDerivative(mesh, tailT, "y");
		if(debug) {
			plotVector(mesh, T_laplace, "T_laplace.dat");
			plotVector(mesh, T_x, "T_x.dat");
			plotVector(mesh, T_y, "T_y.dat");
		}
		Vector q_0 = new SparseVector(dim);
		for(int i=0;i<N-1;i++) {
			q[i] = solveGCM(mesh,i,
					s[i],s[i+1],
					sumLaplaceQ,sumQ_x,sumQ_y,
					//sumLaplaceQ_real,sumQ_x_real,sumQ_y_real,
					T_laplace,T_x,T_y,
					new Vector2Function(phi[i]),
					q_0,phi[i],
					false //flag, if linearize equation q_n
					);
		
			q_0 = q[i];
			
			//sum_{j=1}^{N-1} \Laplace{q}
			sumLaplaceQ.axpy(1.0, Tools.computeLaplace2D(mesh, q[i]));
			sumQ_x.axpy(1.0, Tools.computeDerivative(mesh, q[i], "x"));
			sumQ_y.axpy(1.0, Tools.computeDerivative(mesh, q[i], "y"));//bugfix add 2011-5-7
			if(debug) {
				plotVector(mesh, sumLaplaceQ, "GCM_sumLaplaceQ_"+i+".dat");
				plotVector(mesh, sumQ_x, "GCM_sumQ_x_"+i+".dat");
				plotVector(mesh, sumQ_y, "GCM_sumQ_y_"+i+".dat");
			}
			
			//!!!TEST!!! Real q[i] = phi[i]
//			sumLaplaceQ_real.axpy(1.0,  Tools.computeLaplace2D(mesh, phi[i]));
//			sumQ_x_real.axpy(1.0,  Tools.computeDerivative(mesh, phi[i], "x"));
//			sumQ_y_real.axpy(1.0,  Tools.computeDerivative(mesh, phi[i], "y"));
//			plotVector(mesh, sumLaplaceQ_real, "GCM_sumLaplaceQ_real_"+i+".dat");
//			plotVector(mesh, sumQ_x_real, "GCM_sumQ_x_real_"+i+".dat");
//			plotVector(mesh, sumQ_y_real, "GCM_sumQ_y_real_"+i+".dat");
			
		}
		
		//用最近处光源重构(注意this.delta的位置)：v_tidle -> u ->　a
		Vector v_tidle = new SparseVector(dim);
		for(int i=0;i<N-1;i++) {
			double h = s[i] - s[i+1];
			v_tidle.add(-s[N-1]*s[N-1]*h, q[i]);
			plotVector(mesh, v_tidle, String.format("v_tidle_N_%02d.dat",i));
		}
		v_tidle.add(s[N-1]*s[N-1], tailT);
		plotVector(mesh, v_tidle, "v_tidle_N.dat");
		
		Vector u = new SparseVector(dim);
		for(int i=1;i<=dim;i++) {
			u.set(i, Math.pow(Math.E, v_tidle.get(i)));
		}
		plotVector(mesh, u, "u_N.dat");
		
		Vector a = solveParamInverse(mesh, u);
		plotVector(mesh, a, "alpha_recon_N.dat");
		
		
		//模拟计算出来的phi除了提供边界条件，区域内部的值也有，可以用来重构真实的a(x)
		//*** v_tidle_real:= ln(u) ***
		Vector v_tidle_real = new SparseVector(dim);
		for(int i=0;i<N-1;i++) {
			double h = s[i] - s[i+1];
			v_tidle_real.add(-s[N-1]*s[N-1]*h, phi[i]);
			plotVector(mesh, v_tidle_real, String.format("v_tidle_N_real_%02d.dat",i));
		}
		v_tidle_real.add(s[N-1]*s[N-1], tailT);
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
			Vector phi, //for check only
			boolean linearize
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
//		mapNTF.put(NodeType.Robin, null);
		
		double h = sn_1 - sn; //错误的写成：sn-1 - sn，导致结果异常，以后变量名称一定不要这样写了
		Vector A = getGCMCoef(sn_1,sn);
		if(debug) {
			System.out.println("----------");
			for(int i=1;i<=A.getDim();i++) {
				System.out.print("A"+i+"="+A.get(i)+"  ");
			}
			System.out.println("\n----------");
		}
		//Right hand side of equation q_n
		int dim = sumLaplaceQ.getDim();
		Vector f1 = new SparseVector(dim); //zero vector
		//??? + - ??? A.get(3)*h  //2011-5-7 "-" OK
		f1.add(-A.get(3)*h, sumLaplaceQ);
		f1.add(A.get(3), laplaceT);
		Vector f2x = new SparseVector(dim); //zero vector
		f2x.add(h, sumQ_x);
		f2x.add(-1.0, T_x);
		f2x = FMath.axMuly(A.get(4), f2x, f2x);
		if(debug) plotVector(mesh, f2x, "GCM_f2x_"+n_1+".dat");
		Vector f2y = new SparseVector(dim); //zero vector
		f2y.add(h, sumQ_y);
		f2y.add(-1.0, T_y);
		f2y = FMath.axMuly(A.get(4), f2y, f2y);
		if(debug) plotVector(mesh, f2y, "GCM_f2y_"+n_1+".dat");
		Vector f = FMath.axpy(1.0, f2x, f2y);
		if(debug) plotVector(mesh, f, "GCM_f2xpf2y_"+n_1+".dat");
		f.add(f1);
		if(debug) plotVector(mesh, f, "GCM_"+n_1+"_f.dat");
		
		//b1, b2
		Vector q_x_k_1 = null;
		Vector q_y_k_1 = null;
		
		//------------------------Nonlinear iteration-----------------------
		//Vector qn_1 = q_0;
		Vector qn_1 = new SparseVector(dim);
		Vector qn = null;
		int nMaxIter = 15;
		if(linearize) nMaxIter = 1;
		for(int iter=1;iter<=nMaxIter;iter++) {
			
			//\nabla{q_n} term coefficient (include linearized nonlinear-term) 
			q_x_k_1 = Tools.computeDerivative(mesh, qn_1, "x");
			q_y_k_1 = Tools.computeDerivative(mesh, qn_1, "y");
			if(debug) {
				plotVector(mesh, q_x_k_1, "GCM_q_x_km1_"+n_1+"_"+iter+".dat");
				plotVector(mesh, q_y_k_1, "GCM_q_y_km1_"+n_1+"_"+iter+".dat");
			}
			Vector b1 = new SparseVector(dim); //zero vector
			b1.add(A.get(2)*h, sumQ_x);
			b1.add(-A.get(2), T_x);
			if(!linearize) b1.add(-A.get(1), q_x_k_1);
			Vector b2 = new SparseVector(dim); //zero vector
			b2.add(A.get(2)*h, sumQ_y);
			b2.add(-A.get(2), T_y);
			if(!linearize) b2.add(-A.get(1), q_y_k_1);
			if(debug) {
				plotVector(mesh, b1, "GCM_"+n_1+"_b1_"+iter+".dat");
				plotVector(mesh, b2, "GCM_"+n_1+"_b2_"+iter+".dat");
			}
			//check lhs with rhs
			if(debug && phi !=  null) {
				Vector phi_laplace = Tools.computeLaplace2D(mesh, phi);
				Vector phi_x = Tools.computeDerivative(mesh, phi, "x");
				Vector phi_y = Tools.computeDerivative(mesh, phi, "y");
				plotVector(mesh, phi_laplace, "GCM_"+n_1+"_phi_laplace_"+iter+".dat");
				plotVector(mesh, phi_x, "GCM_"+n_1+"_phi_x_"+iter+".dat");
				plotVector(mesh, phi_y, "GCM_"+n_1+"_phi_y_"+iter+".dat");
				Vector b1r = new SparseVector(dim); //zero vector
				b1r.add(A.get(2)*h, sumQ_x);
				b1r.add(-A.get(2),  T_x);
				b1r.add(-A.get(1),  phi_x);
				Vector b2r = new SparseVector(dim); //zero vector
				b2r.add(A.get(2)*h, sumQ_y);
				b2r.add(-A.get(2),  T_y);
				b2r.add(-A.get(1),  phi_y);
				Vector lhs_x = FMath.axMuly(1.0, phi_x, b1r);
				Vector lhs_y = FMath.axMuly(1.0, phi_y, b2r);
				Vector lhs = FMath.axpy(1.0, lhs_x, lhs_y);
				plotVector(mesh, lhs, "GCM_"+n_1+"_lhsx__lhsy_"+iter+".dat");
				lhs.add(phi_laplace);
				plotVector(mesh, lhs, "GCM_"+n_1+"_lhs_"+iter+".dat");
				//!!!TEST!!! Real \nabla{q_n} term coefficient 
				//b1 = b1r;
				//b2 = b2r;
			}
			
			WeakFormGCM weakForm = new WeakFormGCM();
			weakForm.setF(new Vector2Function(f));
			weakForm.setParam(
					new FC(-1.0),//-1 注意，有负号!!!
					new FC(0.0), //-\eps{q_n}!!! \eps-0.1
					new Vector2Function(b1),
					new Vector2Function(b2)
				);
			
			mesh.clearBorderNodeMark();
			mesh.markBorderNode(mapNTF);

			AssemblerScalar assembler = new AssemblerScalar(mesh, weakForm);
			System.out.println("Begin Assemble...solveGCM");
			assembler.assemble();
			Matrix stiff = assembler.getStiffnessMatrix();
			Vector load = assembler.getLoadVector();
			//Dirichlet condition
			assembler.imposeDirichletCondition(diri);
			System.out.println("Assemble done!");

			SolverJBLAS solver = new SolverJBLAS();
			qn = solver.solveDGESV(stiff, load);
			plotVector(mesh, qn, String.format("GCM_qn_%02d_%02d.dat",n_1,iter));
			
			if(debug) {
				Vector q_x_k = Tools.computeDerivative(mesh, qn, "x");
				Vector q_y_k = Tools.computeDerivative(mesh, qn, "y");
				Vector lhs_x = FMath.axMuly(1.0, q_x_k, b1);
				Vector lhs_y = FMath.axMuly(1.0, q_y_k, b2);
				Vector lhs = FMath.axpy(1.0, lhs_x, lhs_y);
				plotVector(mesh, lhs, "GCM_"+n_1+"_lhsx__lhsy_cmp_"+iter+".dat");
				Vector qn_laplace = Tools.computeLaplace2D(mesh, qn);
				lhs.add(qn_laplace);
				plotVector(mesh, lhs, "GCM_"+n_1+"_lhs_cmp_"+iter+".dat");
				Vector qn_equation = FMath.axpy(-1.0, f, lhs);
				plotVector(mesh, qn_equation, "GCM_"+n_1+"_qn_equation_cmp_"+iter+".dat");
			}
			
			//Stop criterion
			double eps = 0.01;
			System.out.println("GCM Iteration "+iter);
			if(iter > 1) {
				double norm = FMath.axpy(-1.0, qn_1, qn).norm2();
				if(norm < eps) {
					System.out.println("GCM norm < 1.0, STOP");
					break;
				}
			}
			qn_1 = qn;
		}
		return qn;
	}
	
}
