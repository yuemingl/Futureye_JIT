package edu.uta.futureye.application;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;

import edu.uta.futureye.algebra.SparseVectorHashMap;
import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.SparseMatrix;
import edu.uta.futureye.algebra.intf.SparseVector;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.algebra.solver.Solver;
import edu.uta.futureye.algebra.solver.external.SolverJBLAS;
import edu.uta.futureye.core.EdgeLocal;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.core.NodeType;
import edu.uta.futureye.core.intf.AssemblerOld;
import edu.uta.futureye.function.MultiVarFunc;
import edu.uta.futureye.function.FMath;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.basic.FDelta;
import edu.uta.futureye.function.basic.Vector2MathFunc;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.io.MeshReader;
import edu.uta.futureye.lib.assembler.AssemblerScalar;
import edu.uta.futureye.lib.element.FEBilinearRectangleRegular;
import edu.uta.futureye.lib.element.FELinearTriangleOld;
import edu.uta.futureye.lib.weakform.WeakFormL22D;
import edu.uta.futureye.lib.weakform.WeakFormLaplace2D;
import edu.uta.futureye.util.Constant;
import edu.uta.futureye.util.FutureyeException;
import edu.uta.futureye.util.PairDoubleInteger;
import edu.uta.futureye.util.Utils;
import edu.uta.futureye.util.container.ElementList;
import edu.uta.futureye.util.container.NodeList;
import edu.uta.futureye.util.container.ObjList;

public class MouseHead {
	public enum TailType {top,bottom,left,right};
	
	public String outputFolder = "MouseHead";
	public String borderDataPath = null;
	
	//Light source
	public MathFunc delta = null;
	public Variable lightSource = null; //light source position
	
	//Inclusion mu_a
	public MathFunc mu_a = null;
	public double mu_a_bk = 0.1;

	//Inclusion 1/(3*mu_s') = 1.0/30.0 ?
	public MathFunc k = new FC(1.0/50.0);
	
	public double factor = 10000;
	public double baseline = 0.0;
	public double[] factors = null;
	public double[] baselines = null;
	
	//第一种构造tail的方式：分别计算两个反问题，结果相加
	ObjList<Vector> tails = new ObjList<Vector>();
	ObjList<String> tailTypes = new ObjList<String>();
	
	//第二种构造tail的方式：叠加两个外问题的边界延拓，计算一次反问题
	ObjList<Vector> tails2 = new ObjList<Vector>();
	
	public static String gridName = "mouse"; // "mouse","mouse2",..."mouse5"
	
	public static boolean debug = false;
	public static boolean outputMiddleData = false;
	public boolean bSimulate = false;
	public MathFunc mu_aSimulate = null;

	//=true 使用base line测量数据做为背景（通过求解一个外问题）
	//=false 使用计算出来的背景解
	public boolean bUseBaseLineBackground = false;
	
	protected int currentTimeId;
	
	public void setDelta(double x,double y) {
		this.lightSource = new Variable();
		this.lightSource.set("x", x);
		this.lightSource.set("y", y);
		delta = new FDelta(this.lightSource,0.01,2e5);
	}
	
	/**
	 * 给定x坐标位置，构造一个带状的mu_a
	 * 
	 * @param incX
	 * @param incBand
	 * @param maxMu_a
	 * @param bkMu_a
	 */
	public MathFunc makeBandMu_a(double incX, double incBand, double maxMu_a, double bkMu_a) {
		MathFunc rltMu_a = null;
		final double fcx = incX;
		final double fcr = incBand;
		final double fMaxMu_a = maxMu_a;
		final double fBkMu_a = bkMu_a;
		mu_a = new MultiVarFunc("x","y"){
			@Override
			public double apply(Variable v) {
				double dx = v.get("x")-fcx;
				if(Math.sqrt(dx*dx) < fcr) {
					double r = fMaxMu_a*Math.cos((Math.PI/2)*Math.sqrt(dx*dx)/fcr); 
					return r<fBkMu_a?fBkMu_a:r;
				}
				else
					return fBkMu_a;
			}
		};
		return rltMu_a;
	}

	/**
	 * 构造一个圆形mu_a
	 * 
	 * @param incX
	 * @param incY
	 * @param incR
	 * @param maxMu_a
	 * @return
	 */
	public MathFunc makeMu_a(double incX, double incY, double incR, double maxMu_a, double bkMu_a) {
		MathFunc rltMu_a = null;
		final double fcx = incX;
		final double fcy = incY;
		final double fcr = incR;
		final double fMaxMu_a = maxMu_a;
		final double fBkMu_a = bkMu_a;
		rltMu_a = new MultiVarFunc("x","y"){
			@Override
			public double apply(Variable v) {
				double dx = v.get("x")-fcx;
				double dy = v.get("y")-fcy;
				if(Math.sqrt(dx*dx+dy*dy) < fcr) {
					double r = fMaxMu_a*Math.cos((Math.PI/2)*Math.sqrt(dx*dx+dy*dy)/fcr); 
					return r<fBkMu_a?fBkMu_a:r;
				}
				else
					return fBkMu_a;
			}
		};
		return rltMu_a;
	}
	
	public Vector solveForwardNeumann(Mesh mesh) {
		//Mark border type
		HashMap<NodeType, MathFunc> mapNTF = new HashMap<NodeType, MathFunc>();
		mapNTF.clear();
		mapNTF.put(NodeType.Robin, null);
		mesh.clearBorderNodeMark();
		mesh.markBorderNode(mapNTF);
		
		WeakFormLaplace2D weakForm = new WeakFormLaplace2D();
		
		//Right hand side
		weakForm.setF(this.delta);
		
		// *** u + u_n = 0, on boundary ***
		weakForm.setParam(
				this.k, this.mu_a, FC.C0, this.k //d==k,q=0 (即：u_n + u =0)
			);
		
		AssemblerOld assembler = new AssemblerScalar(mesh, weakForm);
		System.out.println("Begin Assemble...solveForwardNeumann");
		assembler.assemble();
		SparseMatrix stiff = assembler.getStiffnessMatrix();
		SparseVector load = assembler.getLoadVector();
		//assembler.imposeDirichletCondition(new FC(0.0));
		System.out.println("Assemble done!");

		Solver solver = new Solver();
		Vector u = solver.solveCGS(stiff, load);
		return u;
	}	
	
	public Vector solveForwardDirichlet(Mesh mesh, MathFunc diri) {
		WeakFormLaplace2D weakForm = new WeakFormLaplace2D();
		
		//Right hand side
		weakForm.setF(this.delta);

		weakForm.setParam(
				this.k, this.mu_a, FC.C0, this.k //d==k,q=0 (即：u_n + u =0)
			);
		
		//!!!AssemblerScalar与AssemblerScalarFast有很大区别，问题？？？ 20110921
		//Assembler assembler = new AssemblerScalarFast(mesh, weakForm);
		AssemblerOld assembler = new AssemblerScalar(mesh, weakForm);
		System.out.println("Begin Assemble...solveForwardDirichlet");
		assembler.assemble();
		SparseMatrix stiff = assembler.getStiffnessMatrix();
		SparseVector load = assembler.getLoadVector();
		//Dirichlet condition
		assembler.imposeDirichletCondition(diri);
		
		System.out.println("Assemble done!");

		Solver solver = new Solver();
		Vector u = solver.solveCGS(stiff, load);
		return u;
	}	
	
	//public Vector solveParamInverse(Mesh mesh, Vector U, final TailType type) {
	public Vector solveParamInverse(Mesh mesh, Vector U) {
		HashMap<NodeType, MathFunc> mapNTF2 = new HashMap<NodeType, MathFunc>();
		
		mapNTF2.put(NodeType.Dirichlet, new MultiVarFunc("x","y") {
			@Override
			public double apply(Variable v) {
				//应该全是Dirichlet条件，否则在非Dirichlet条件处会产生大变化
//				double x = v.get("x");
//				double y = v.get("y");
//				if(type == TailType.left) {
//					if(Math.abs(x-(-1.6))<Constant.meshEps) {
//						return 0;
//					} else 
//						return 1;
//				}
//				if(type == TailType.right) {
//					if(Math.abs(x-(1.6))<Constant.meshEps) {
//						return 0;
//					} else 
//						return 1;
//				}
//				if(type == TailType.top) {
//					if(Math.abs(y-(0.9))<Constant.meshEps) {
//						return 0;
//					} else 
//						return 1;
//				}
//				if(type == TailType.bottom) {
//					if(Math.abs(y-(-0.5))<Constant.meshEps) {
//						return 0;
//					} else 
//						return 1;
//				}
				return 1;
			}			
		});
		mesh.clearBorderNodeMark();
		mesh.markBorderNode(mapNTF2);
		
		//Weak form
		WeakFormL22D weakFormL2 = new WeakFormL22D();
		//Right hand side
		weakFormL2.setF(this.delta);
		//Parameters
		weakFormL2.setParam(
				this.k, new Vector2MathFunc(U)
			);
		
		AssemblerOld assembler = new AssemblerScalar(mesh, weakFormL2);
		System.out.println("Begin Assemble...solveParamInverse");
		assembler.assemble();
		Matrix stiff = assembler.getStiffnessMatrix();
		Vector load = assembler.getLoadVector();
		assembler.imposeDirichletCondition(new FC(this.mu_a_bk));
		System.out.println("Assemble done!");
		
		SolverJBLAS solver = new SolverJBLAS();
		Vector u = solver.solveDGESV(stiff, load);
		return u;
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
		Vector rlt = new SparseVectorHashMap(dimTo);
		for(int i=1;i<=nodeTo.size();i++) {
			Node node = meshFrom.findNode(nodeTo.at(i));
			if(node != null) {
				rlt.set(nodeTo.at(i).globalIndex, u.get(node.globalIndex));
			}
		}
		return rlt;
	}
	
	/**
	 * 提取网格Omega1（外问题）内边界点的参数坐标，提供给matlab程序的
	 * “插值结点的参数坐标(0->r1->r2->r3->r4(0))”
	 * 
	 *  S1            S6
	 * (r4)0---->----r1
	 *     /          \
	 *    /            \
	 *   ┐              ┘
	 *  /                \
	 * r3-------<--------r2
	 * S4                S9
	 * 
	 * 在标准输出打印结点编号、参数坐标对列表，以参数坐标值顺序排序
	 */	
	public void extractBorderParamData(Mesh omega) {
		ObjList<PairDoubleInteger> rs = new ObjList<PairDoubleInteger>();
		NodeList nodes = omega.getNodeList();
		Node S1 = new Node(0,-0.9,0.8);
		Node S6 = new Node(0,0.9,0.8);
		Node S9 = new Node(0,1.5,-0.4);
		Node S4 = new Node(0,-1.5,-0.4);
		double dis1 = Utils.computeLength(S1,S6);
		double dis2 = dis1 + Utils.computeLength(S6,S9);
		double dis3 = dis2 + Utils.computeLength(S9,S4);
		//double dis4 = dis3 + Utils.computeLength(S4,S1);
		
		for(int i=1;i<=nodes.size();i++) {
			Node node = nodes.at(i);
			if(node.isInnerNode())
				continue;
			//除去小矩形（外问题在Omega1上求解时）
			if(Math.abs(node.coord(1)-(-1.6))<Constant.meshEps ||
					Math.abs(node.coord(1)-1.6)<Constant.meshEps ||
					Math.abs(node.coord(2)-(-0.5))<Constant.meshEps ||
					Math.abs(node.coord(2)-0.9)<Constant.meshEps)
				continue;
			//除去大矩形（外问题在Omega11上求解）
			if(Math.abs(node.coord(1)-(-2.1))<Constant.meshEps ||
					Math.abs(node.coord(1)-2.1)<Constant.meshEps ||
					Math.abs(node.coord(2)-(-1.0))<Constant.meshEps ||
					Math.abs(node.coord(2)-1.4)<Constant.meshEps)
				continue;
			
			PairDoubleInteger di = new PairDoubleInteger();
			di.i = node.getIndex();
			
			if(Utils.isPointOnLineSegment(S1, S6, node)) {//top
				di.d = Utils.computeLength(S1,node);
			} else if (Utils.isPointOnLineSegment(S6, S9, node)) {//right
				di.d = dis1 + Utils.computeLength(S6,node);
			} else if(Utils.isPointOnLineSegment(S9, S4, node)) {//bottom
				di.d = dis2 + Utils.computeLength(S9,node);
			} else if(Utils.isPointOnLineSegment(S4, S1, node)) {//left
				di.d = dis3 + Utils.computeLength(S4,node);
			} else {
				throw new FutureyeException("");
			}
			rs.add(di);
		}
		
		Collections.sort(rs.toList(), new Comparator<PairDoubleInteger>() {
			@Override
			public int compare(PairDoubleInteger o1, PairDoubleInteger o2) {
				if(o1.d>o2.d) return 1;
				else return -1;
			}
			
		});
		for(int i=1;i<=rs.size();i++) {
			System.out.println(rs.at(i).i+" "+rs.at(i).d);
		}
	}
	
	/**
	 * 读取matlab插值后的数据：(globalIndex, value)对
	 * @param srcId
	 * @return
	 */
	public ObjList<PairDoubleInteger> readInterpData(int srcId, int timeId) {
		ObjList<PairDoubleInteger> vs = new ObjList<PairDoubleInteger>();
		FileInputStream in;
		try {
			in = new FileInputStream(
					String.format("./"+outputFolder+"/"+borderDataPath+"/InterpOut%d_%d.txt", 
							srcId,timeId));

			InputStreamReader reader = new InputStreamReader(in);
			BufferedReader br = new BufferedReader(reader);
	
			String str = null;
			while((str = br.readLine()) != null){
				if(debug)
					System.out.println(str);
				if(str.startsWith("#")) continue;
				String[] line = str.split("(\\s)+");
				PairDoubleInteger pair = new PairDoubleInteger();
				pair.i = Integer.parseInt(line[0]);
				if(factors != null)
					pair.d = Double .parseDouble(line[2])*factors[srcId] + baselines[srcId];
				else
					pair.d = Double .parseDouble(line[2])*factor + baseline;
				vs.add(pair);
			}
			br.close();
			in.close();
			return vs;
		} catch (Exception e) {
			e.printStackTrace();
		}
		return null;
	}
	
	
	/**
	 * 获取区域的边界结点列表（与Omega2边界一致的任何区域）
	 * Omega2:  [-1.6,1.6]*[-0.5,0.9]
	 * 
	 * @param type
	 * @param omega2
	 * @return
	 */
	public ObjList<Node> getBorderNodes_Omega2(TailType type, Mesh omega2) {
		NodeList nodes2 = omega2.getNodeList(); 
		ObjList<Node> borderNodes = new ObjList<Node>();
		double x,y;
		if(type == TailType.left) {
			for(int i=1;i<=nodes2.size();i++) {
				if( //y=[-0.5,0.9]
					nodes2.at(i).coord(2)-(-0.5)>-Constant.meshEps && 0.9-nodes2.at(i).coord(2)>-Constant.meshEps &&
					//x=-1.6
					Math.abs(nodes2.at(i).coord(1)-(-1.6))<Constant.meshEps
				) {
					borderNodes.add(nodes2.at(i));
					if(debug)
						System.out.println(nodes2.at(i));
				}
			}
		} else if(type == TailType.bottom) {
			for(int i=1;i<=nodes2.size();i++) {
				x = nodes2.at(i).coord(1);
				y = nodes2.at(i).coord(2);
				if( //x:[-1.6,1.6]
					x>-1.6-Constant.meshEps && x<1.6+Constant.meshEps &&
					//y:-0.5
					Math.abs(y-(-0.5))<Constant.meshEps
				) {
					borderNodes.add(nodes2.at(i));
					if(debug)
						System.out.println(nodes2.at(i));
				}
			}
		}  else if(type == TailType.right) {
			for(int i=1;i<=nodes2.size();i++) {
				if( //y=[-0.5,0.9]
					nodes2.at(i).coord(2)-(-0.5)>-Constant.meshEps && 0.9-nodes2.at(i).coord(2)>-Constant.meshEps &&
					//x=1.6
					Math.abs(nodes2.at(i).coord(1)-1.6)<Constant.meshEps
					) {
						borderNodes.add(nodes2.at(i));
						if(debug)
							System.out.println(nodes2.at(i));
					}
			}
		}  else if(type == TailType.top) {
			for(int i=1;i<=nodes2.size();i++) {
				if( //x:[-1.6,1.6]
					nodes2.at(i).coord(1)-(-1.6)>-Constant.meshEps && 1.6-nodes2.at(i).coord(1)>-Constant.meshEps &&
					//y:0.9
					Math.abs(nodes2.at(i).coord(2)-(0.9))<Constant.meshEps
				) {
					borderNodes.add(nodes2.at(i));
					if(debug)
						System.out.println(nodes2.at(i));
				}
			}
		} else {
			throw new FutureyeException(""+type);
		}
		return borderNodes;
	}
	
	/**
	 * 计算光源在右边的tail=>左tail
	 */	
	public Vector computeTailLeft(Mesh omega2,Vector u2_bk,
			Mesh omega1,Vector u1) {
		Vector tail = new SparseVectorHashMap(u2_bk.getDim());
		NodeList nodes2 = omega2.getNodeList(); 
		ObjList<Node> u2LeftNodes = this.getBorderNodes_Omega2(TailType.left, omega2);
		ObjList<Node> u1LeftNodes = this.getBorderNodes_Omega2(TailType.left, omega1);
		
		//on left, tail=u1-u0
		for(int i=1;i<=u2LeftNodes.size();i++) {
			Node u2BottomNode = u2LeftNodes.at(i);
			for(int j=1;j<=u1LeftNodes.size();j++) {
				Node u1BottomNode = u1LeftNodes.at(j);
				if(u2BottomNode.coordEquals(u1BottomNode)) {
					tail.set(u2BottomNode.globalIndex, 
							u1.get(u1BottomNode.globalIndex)-u2_bk.get(u2BottomNode.globalIndex));
				}
			}
		}
		Tools.plotVector(omega2, outputFolder, "left_tail0.dat", tail);
		
		//extends from left to omega2
		for(int i=1;i<=nodes2.size();i++) {
			Node node  = nodes2.at(i);
			
			if(Math.abs(node.coord(1)-(-1.6))<Constant.meshEps) 
				continue;
			
			boolean find = false;
			for(int j=1;j<=u2LeftNodes.size();j++) {
				if(find) break;
				Node u2LeftNode = u2LeftNodes.at(j);
				ElementList eles = u2LeftNode.belongToElements;
				for(int k=1;k<=eles.size();k++) {
					if(find) break;
					Element e = eles.at(k);
					ObjList<EdgeLocal> edges = e.edges();
					for(int ie=1;ie<=edges.size();ie++) {
						if(find) break;
						EdgeLocal edge = edges.at(ie);
						Node newNode = new Node(0,u2LeftNode.coord(1),node.coord(2));
						if(Utils.isPointOnLineSegment(edge.beginNode(), edge.endNode(),
								newNode)) {
							double val = Utils.linearInterpolate(edge.beginNode(), edge.endNode(),
									newNode, 
									tail.get(edge.beginNode().globalIndex), 
									tail.get(edge.endNode().globalIndex));
							//光源
							double lx = lightSource.get("x");
							double ly = lightSource.get("y");
							//内部点
							double x = node.coord(1);
							double y = node.coord(2);
							
							double MaxS = Math.sqrt((lx-(-1.6))*(lx-(-1.6)));
							double LS = Math.sqrt((x-lx)*(x-lx)+(y-ly)*(y-ly));
							if(LS>MaxS) LS=MaxS;
							double Inc = Math.exp(3.3*(MaxS-LS)/MaxS)-1.0;
							
							double SS = val*(1+Inc);
							//tail.set(node.globalIndex, SS);
							tail.set(node.globalIndex, val);
							find=true;
						}
					}
				}
			}
		}
		Tools.plotVector(omega2, outputFolder, String.format("left_tail1_%03d.dat",currentTimeId), tail);
		tails2.add(tail.copy());
//		tail = Utils.gaussSmooth(omega2, tail, 2, 0.5);
//		tail = Utils.gaussSmooth(omega2, tail, 2, 0.4);
//		tail = Utils.gaussSmooth(omega2, tail, 2, 0.3);
//		tail = Utils.gaussSmooth(omega2, tail, 2, 0.3);
//		tail = Utils.gaussSmooth(omega2, tail, 2, 0.3);
//		tail = Utils.gaussSmooth(omega2, tail, 2, 0.3);
		
		Tools.plotVector(omega2, outputFolder, "left_tail1_smooth.dat", tail);
		
		for(int i=1;i<=nodes2.size();i++) {
			Node node  = nodes2.at(i);
			tail.add(node.globalIndex,u2_bk.get(node.globalIndex));
		}
		Tools.plotVector(omega2, outputFolder, "left_tail.dat", tail);
		
		
		return tail;
	}
	
	public Vector computeTailRight(Mesh omega2,Vector u2_bk,
			Mesh omega1,Vector u1) {
		Vector tail = new SparseVectorHashMap(u2_bk.getDim());
		NodeList nodes2 = omega2.getNodeList(); 
		ObjList<Node> u2RightNodes = this.getBorderNodes_Omega2(TailType.right, omega2);
		ObjList<Node> u1RightNodes = this.getBorderNodes_Omega2(TailType.right, omega1);
		
		//on left tail=u1-u0
		for(int i=1;i<=u2RightNodes.size();i++) {
			Node u2BottomNode = u2RightNodes.at(i);
			for(int j=1;j<=u1RightNodes.size();j++) {
				Node u1BottomNode = u1RightNodes.at(j);
				if(u2BottomNode.coordEquals(u1BottomNode)) {
					tail.set(u2BottomNode.globalIndex, 
							u1.get(u1BottomNode.globalIndex)-u2_bk.get(u2BottomNode.globalIndex));
				}
			}
		}
		Tools.plotVector(omega2, outputFolder, "right_tail0.dat", tail);
		
		//extends from right to omega2
		for(int i=1;i<=nodes2.size();i++) {
			Node node  = nodes2.at(i);
			
			if(Math.abs(node.coord(1)-1.6)<Constant.meshEps) 
				continue;
			
			boolean find = false;
			for(int j=1;j<=u2RightNodes.size();j++) {
				if(find) break;
				Node u2RightNode = u2RightNodes.at(j);
				ElementList eles = u2RightNode.belongToElements;
				for(int k=1;k<=eles.size();k++) {
					if(find) break;
					Element e = eles.at(k);
					ObjList<EdgeLocal> edges = e.edges();
					for(int ie=1;ie<=edges.size();ie++) {
						if(find) break;
						EdgeLocal edge = edges.at(ie);
						Node newNode = new Node(0,u2RightNode.coord(1),node.coord(2));
						if(Utils.isPointOnLineSegment(edge.beginNode(), edge.endNode(),
								newNode)) {
							double val = Utils.linearInterpolate(edge.beginNode(), edge.endNode(),
									newNode, 
									tail.get(edge.beginNode().globalIndex), 
									tail.get(edge.endNode().globalIndex));
							//光源
							double lx = lightSource.get("x");
							double ly = lightSource.get("y");
							//内部点
							double x = node.coord(1);
							double y = node.coord(2);
							
							double MaxS = Math.sqrt((lx-1.6)*(lx-1.6));
							double LS = Math.sqrt((x-lx)*(x-lx)+(y-ly)*(y-ly));
							if(LS>MaxS) LS=MaxS;
							double Inc = Math.exp(3.3*(MaxS-LS)/MaxS)-1.0;
							
							double SS = val*(1+Inc);
							//tail.set(node.globalIndex, SS);
							tail.set(node.globalIndex, val);
							find=true;
						}
					}
				}
			}
		}
		Tools.plotVector(omega2, outputFolder, String.format("right_tail1_%03d.dat",currentTimeId), tail);
		tails2.add(tail.copy());
		
		
		for(int i=1;i<=nodes2.size();i++) {
			Node node  = nodes2.at(i);
			tail.add(node.globalIndex,u2_bk.get(node.globalIndex));
		}
		Tools.plotVector(omega2, outputFolder, "right_tail.dat", tail);
		
		return tail;
	}
	
	/**
	 * 计算光源在上边的tail=>下tail
	 */
	public Vector computeTailBottom(Mesh omega2,Vector u2_bk,
			Mesh omega1,Vector u1) {
		Vector tail = new SparseVectorHashMap(u2_bk.getDim());
		NodeList nodes2 = omega2.getNodeList(); 
		ObjList<Node> u2BottomNodes = this.getBorderNodes_Omega2(TailType.bottom, omega2);
		ObjList<Node> u1BottomNodes = this.getBorderNodes_Omega2(TailType.bottom, omega1);
		
		//on bottom tail=u1-u0
		for(int i=1;i<=u2BottomNodes.size();i++) {
			Node u2BottomNode = u2BottomNodes.at(i);
			for(int j=1;j<=u1BottomNodes.size();j++) {
				Node u1BottomNode = u1BottomNodes.at(j);
				if(u2BottomNode.coordEquals(u1BottomNode)) {
					tail.set(u2BottomNode.globalIndex, 
							u1.get(u1BottomNode.globalIndex)-u2_bk.get(u2BottomNode.globalIndex));
				}
			}
		}
		Tools.plotVector(omega2, outputFolder, "bottom_tail0.dat", tail);
		
		//extends from bottom to omega2
		for(int i=1;i<=nodes2.size();i++) {
			Node node  = nodes2.at(i);
			
			if(Math.abs(node.coord(2)-(-0.5))<Constant.meshEps) 
				continue;
			
			boolean find = false;
			for(int j=1;j<=u2BottomNodes.size();j++) {
				if(find) break;
				Node u2BottomNode = u2BottomNodes.at(j);
				ElementList eles = u2BottomNode.belongToElements;
				for(int k=1;k<=eles.size();k++) {
					if(find) break;
					Element e = eles.at(k);
					ObjList<EdgeLocal> edges = e.edges();
					for(int ie=1;ie<=edges.size();ie++) {
						if(find) break;
						EdgeLocal edge = edges.at(ie);
						Node newNode = new Node(0,node.coord(1),u2BottomNode.coord(2));
						if(Utils.isPointOnLineSegment(edge.beginNode(), edge.endNode(),
								newNode)) {
							double val = Utils.linearInterpolate(edge.beginNode(), edge.endNode(),
									newNode, 
									tail.get(edge.beginNode().globalIndex), 
									tail.get(edge.endNode().globalIndex));
							//光源
							double lx = lightSource.get("x");
							double ly = lightSource.get("y");
							//内部点
							double x = node.coord(1);
							double y = node.coord(2);
							
							double MaxS = Math.sqrt((ly-(-0.5))*(ly-(-0.5)));
							double LS = Math.sqrt((x-lx)*(x-lx)+(y-ly)*(y-ly));
							if(LS>MaxS) LS=MaxS;
							//double Inc = Math.exp(BS/(LS+BS))-1;
							double Inc = Math.exp(3.3*(MaxS-LS)/MaxS)-1.0;
							
							double SS = val*(1+Inc);
							//tail.set(node.globalIndex, SS);
							tail.set(node.globalIndex, val);
							find=true;
						}
					}
				}
			}
		}
		Tools.plotVector(omega2, outputFolder, String.format("bottom_tail1_%03d.dat",currentTimeId), tail);
		tails2.add(tail.copy());
//		tail = Utils.gaussSmooth(omega2, tail, 2, 0.5);
//		tail = Utils.gaussSmooth(omega2, tail, 2, 0.4);
//		tail = Utils.gaussSmooth(omega2, tail, 2, 0.3);
//		tail = Utils.gaussSmooth(omega2, tail, 2, 0.3);
//		tail = Utils.gaussSmooth(omega2, tail, 2, 0.3);
//		tail = Utils.gaussSmooth(omega2, tail, 2, 0.3);
		
		
		for(int i=1;i<=nodes2.size();i++) {
			Node node  = nodes2.at(i);
			tail.add(node.globalIndex,u2_bk.get(node.globalIndex));
		}
		Tools.plotVector(omega2, outputFolder, "bottom_tail.dat", tail);
		
		return tail;
	}
	
	/**
	 * 计算光源在上边的tail=>下tail
	 */
	public Vector computeTailTop(Mesh omega2,Vector u2_bk,
			Mesh omega1,Vector u1) {
		Vector tail = new SparseVectorHashMap(u2_bk.getDim());
		NodeList nodes2 = omega2.getNodeList(); 
		ObjList<Node> u2TopNodes = this.getBorderNodes_Omega2(TailType.top, omega2);
		ObjList<Node> u1TopNodes = this.getBorderNodes_Omega2(TailType.top, omega1);
		
		//on top tail=u1-u0
		for(int i=1;i<=u2TopNodes.size();i++) {
			Node u2BottomNode = u2TopNodes.at(i);
			for(int j=1;j<=u1TopNodes.size();j++) {
				Node u1BottomNode = u1TopNodes.at(j);
				if(u2BottomNode.coordEquals(u1BottomNode)) {
					tail.set(u2BottomNode.globalIndex, 
							u1.get(u1BottomNode.globalIndex)-u2_bk.get(u2BottomNode.globalIndex));
				}
			}
		}
		Tools.plotVector(omega2, outputFolder, "top_tail0.dat", tail);
		
		//extends from top to omega2
		for(int i=1;i<=nodes2.size();i++) {
			Node node  = nodes2.at(i);
			
			if(Math.abs(node.coord(2)-1.0)<Constant.meshEps) 
				continue;
			
			boolean find = false;
			for(int j=1;j<=u2TopNodes.size();j++) {
				if(find) break;
				Node u2TopNode = u2TopNodes.at(j);
				ElementList eles = u2TopNode.belongToElements;
				for(int k=1;k<=eles.size();k++) {
					if(find) break;
					Element e = eles.at(k);
					ObjList<EdgeLocal> edges = e.edges();
					for(int ie=1;ie<=edges.size();ie++) {
						if(find) break;
						EdgeLocal edge = edges.at(ie);
						Node newNode = new Node(0,node.coord(1),u2TopNode.coord(2));
						if(Utils.isPointOnLineSegment(edge.beginNode(), edge.endNode(),
								newNode)) {
							double val = Utils.linearInterpolate(edge.beginNode(), edge.endNode(),
									newNode, 
									tail.get(edge.beginNode().globalIndex), 
									tail.get(edge.endNode().globalIndex));
							//光源
							double lx = lightSource.get("x");
							double ly = lightSource.get("y");
							//内部点
							double x = node.coord(1);
							double y = node.coord(2);
							
							double MaxS = Math.sqrt((ly-1.0)*(ly-1.0));
							double LS = Math.sqrt((x-lx)*(x-lx)+(y-ly)*(y-ly));
							if(LS>MaxS) LS=MaxS;
							//double Inc = Math.exp(BS/(LS+BS))-1;
							double Inc = Math.exp(3.3*(MaxS-LS)/MaxS)-1.0;
							
							double SS = val*(1+Inc);
							//tail.set(node.globalIndex, SS);
							tail.set(node.globalIndex, val);
							find=true;
						}
					}
				}
			}
		}
		Tools.plotVector(omega2, outputFolder, String.format("top_tail1_%03d.dat",currentTimeId), tail);
		tails2.add(tail.copy());
		
		for(int i=1;i<=nodes2.size();i++) {
			Node node  = nodes2.at(i);
			tail.add(node.globalIndex,u2_bk.get(node.globalIndex));
		}
		Tools.plotVector(omega2, outputFolder, "top_tail.dat", tail);
		
		return tail;
	}

	
	public void extendAlphaFromBottom(Mesh omega2,Vector alpha) {
		NodeList nodes2 = omega2.getNodeList(); 
		ObjList<Node> u2BottomNodes = this.getBorderNodes_Omega2(TailType.bottom, omega2);
		//extends from bottom to omega2
		for(int i=1;i<=nodes2.size();i++) {
			Node node  = nodes2.at(i);
			
			if(Math.abs(node.coord(2)-(-0.5))<Constant.meshEps) 
				continue;
			
			boolean find = false;
			for(int j=1;j<=u2BottomNodes.size();j++) {
				if(find) break;
				Node u2BottomNode = u2BottomNodes.at(j);
				ElementList eles = u2BottomNode.belongToElements;
				for(int k=1;k<=eles.size();k++) {
					if(find) break;
					Element e = eles.at(k);
					ObjList<EdgeLocal> edges = e.edges();
					for(int ie=1;ie<=edges.size();ie++) {
						if(find) break;
						EdgeLocal edge = edges.at(ie);
						Node newNode = new Node(0,node.coord(1),u2BottomNode.coord(2));
						if(Utils.isPointOnLineSegment(edge.beginNode(), edge.endNode(),
								newNode)) {
							double val = Utils.linearInterpolate(edge.beginNode(), edge.endNode(),
									newNode, 
									alpha.get(edge.beginNode().globalIndex), 
									alpha.get(edge.endNode().globalIndex));
							alpha.set(node.globalIndex, val);
							find=true;
						}
					}
				}
			}
		}			
	}
	
	public void extendAlphaFromLeft(Mesh omega2,Vector alpha) {
		NodeList nodes2 = omega2.getNodeList(); 
		ObjList<Node> u2LeftNodes = this.getBorderNodes_Omega2(TailType.left, omega2);
		
		//extends from left to omega2
		for(int i=1;i<=nodes2.size();i++) {
			Node node  = nodes2.at(i);
			
			if(Math.abs(node.coord(1)-(-1.6))<Constant.meshEps) 
				continue;
			
			boolean find = false;
			for(int j=1;j<=u2LeftNodes.size();j++) {
				if(find) break;
				Node u2LeftNode = u2LeftNodes.at(j);
				ElementList eles = u2LeftNode.belongToElements;
				for(int k=1;k<=eles.size();k++) {
					if(find) break;
					Element e = eles.at(k);
					ObjList<EdgeLocal> edges = e.edges();
					for(int ie=1;ie<=edges.size();ie++) {
						if(find) break;
						EdgeLocal edge = edges.at(ie);
						Node newNode = new Node(0,u2LeftNode.coord(1),node.coord(2));
						if(Utils.isPointOnLineSegment(edge.beginNode(), edge.endNode(),
								newNode)) {
							double val = Utils.linearInterpolate(edge.beginNode(), edge.endNode(),
									newNode, 
									alpha.get(edge.beginNode().globalIndex), 
									alpha.get(edge.endNode().globalIndex));
							alpha.set(node.globalIndex, val);
							find=true;
						}
					}
				}
			}
		}
			
	}
	
	public void extendAlphaFromRight(Mesh omega2,Vector alpha) {
		NodeList nodes2 = omega2.getNodeList(); 
		ObjList<Node> u2RightNodes = this.getBorderNodes_Omega2(TailType.right, omega2);
		
		//extends from right to omega2
		for(int i=1;i<=nodes2.size();i++) {
			Node node  = nodes2.at(i);
			
			if(Math.abs(node.coord(1)-1.6)<Constant.meshEps) 
				continue;
			
			boolean find = false;
			for(int j=1;j<=u2RightNodes.size();j++) {
				if(find) break;
				Node u2RightNode = u2RightNodes.at(j);
				ElementList eles = u2RightNode.belongToElements;
				for(int k=1;k<=eles.size();k++) {
					if(find) break;
					Element e = eles.at(k);
					ObjList<EdgeLocal> edges = e.edges();
					for(int ie=1;ie<=edges.size();ie++) {
						if(find) break;
						EdgeLocal edge = edges.at(ie);
						Node newNode = new Node(0,u2RightNode.coord(1),node.coord(2));
						if(Utils.isPointOnLineSegment(edge.beginNode(), edge.endNode(),
								newNode)) {
							double val = Utils.linearInterpolate(edge.beginNode(), edge.endNode(),
									newNode, 
									alpha.get(edge.beginNode().globalIndex), 
									alpha.get(edge.endNode().globalIndex));
							alpha.set(node.globalIndex, val);
							find=true;
						}
					}
				}
			}
		}
	}
	
	public void extendAlphaFromTop(Mesh omega2,Vector alpha) {
		NodeList nodes2 = omega2.getNodeList(); 
		ObjList<Node> u2TopNodes = this.getBorderNodes_Omega2(TailType.top, omega2);
		
		//extends from top to omega2
		for(int i=1;i<=nodes2.size();i++) {
			Node node  = nodes2.at(i);
			
			if(Math.abs(node.coord(2)-1.0)<Constant.meshEps) 
				continue;
			
			boolean find = false;
			for(int j=1;j<=u2TopNodes.size();j++) {
				if(find) break;
				Node u2TopNode = u2TopNodes.at(j);
				ElementList eles = u2TopNode.belongToElements;
				for(int k=1;k<=eles.size();k++) {
					if(find) break;
					Element e = eles.at(k);
					ObjList<EdgeLocal> edges = e.edges();
					for(int ie=1;ie<=edges.size();ie++) {
						if(find) break;
						EdgeLocal edge = edges.at(ie);
						Node newNode = new Node(0,node.coord(1),u2TopNode.coord(2));
						if(Utils.isPointOnLineSegment(edge.beginNode(), edge.endNode(),
								newNode)) {
							double val = Utils.linearInterpolate(edge.beginNode(), edge.endNode(),
									newNode, 
									alpha.get(edge.beginNode().globalIndex), 
									alpha.get(edge.endNode().globalIndex));
							alpha.set(node.globalIndex, val);
							find=true;
						}
					}
				}
			}
		}
	}

	/**
	 * 根据不同的tail标记不同的边界类型
	 * @param mesh
	 * @param fTailType
	 */
	public void markExteriorBorder(Mesh mesh, final TailType fTailType) {
		//外问题 Solve exterior problem
		HashMap<NodeType, MathFunc> mapNTF = new HashMap<NodeType, MathFunc>();
		mapNTF.put(NodeType.Robin, new MultiVarFunc("x","y") {
			@Override
			public double apply(Variable v) {
				double x = v.get("x");
				double y = v.get("y");
				//Omega1
//				if(Math.abs(x-1.6)<Constant.meshEps ||
//						Math.abs(x+1.6)<Constant.meshEps ||
//						Math.abs(y-0.9)<Constant.meshEps ||
//						Math.abs(y+0.5)<Constant.meshEps)
				//Omega11 外问题区域的外边界
				if(Math.abs(x-2.1)<Constant.meshEps ||
						Math.abs(x+2.1)<Constant.meshEps ||
						Math.abs(y-1.4)<Constant.meshEps ||
						Math.abs(y+1.0)<Constant.meshEps)
					return 1;
				//bottom tail: 外问题上边界不要Dirichlet条件，改为Robin
				if(fTailType==TailType.bottom && Math.abs(y-0.8)<Constant.meshEps)
					return 1;
				//top tail: 外问题右边界不要Dirichlet条件，改为Robin
				if(fTailType==TailType.top && Math.abs(y-(-0.4))<Constant.meshEps)
					return 1;
				//最后处理left tail和right tail，两条边界是斜线，只能指定一个范围：
				//left tail: 外问题右边界不要Dirichlet条件，改为Robin
				if(fTailType==TailType.left && Math.abs(x-1.6)<Constant.meshEps) //应该改为：x>0.9-Constant.meshEps
					return 1;
				//right tail: 外问题右边界不要Dirichlet条件，改为Robin
				if(fTailType==TailType.right && Math.abs(x-(-1.6))<Constant.meshEps) //应该改为：x<-0.9+Constant.meshEps
					return 1;
				
				return 0;
			}
		});
		
		mapNTF.put(NodeType.Dirichlet, null);
		mesh.clearBorderNodeMark();
		mesh.markBorderNode(mapNTF);		
	}
	
	public boolean isOnTopBorderOmega(Node node) {
		double x = node.coord(1);
		double y = node.coord(2);
		if(Math.abs(y-0.8)<Constant.meshEps &&
				x>=-0.9-Constant.meshEps && x<=0.9+Constant.meshEps)
			return true;
		else
			return false;
	}
	
	public boolean isOnBottomBorderOmega(Node node) {
		double x = node.coord(1);
		double y = node.coord(2);
		if(Math.abs(y-(-0.4))<Constant.meshEps &&
				x>=-1.5-Constant.meshEps && x<=1.5+Constant.meshEps)
			return true;
		else
			return false;
	}
	
	public boolean isOnLeftBorderOmega(Node node) {
		Node S25 = new Node(0,-0.9,0.8);
		Node S28 = new Node(0,-1.5,-0.4);
		return Utils.isPointOnLineSegment(S25, S28, node);
	}
	
	public boolean isOnRightBorderOmega(Node node) {
		Node S1 = new Node(0,0.9,0.8);
		Node S4 = new Node(0,1.5,-0.4);
		return Utils.isPointOnLineSegment(S1, S4, node);
	}	
	
	public void run(int srcLightId, int timeId, TailType tailType) {
		currentTimeId = timeId;
		
		//Read a triangle mesh from an input file
		//老鼠头，梯形区域
		MeshReader reader = new MeshReader(gridName+"_omega.grd");
		//包围老鼠头的矩形规则区域
		MeshReader reader2 = new MeshReader(gridName+"_omega2.grd");
		//大区域背景解
		MeshReader reader0 = new MeshReader(gridName+"_omega00.grd");
		//外问题，大区域
		MeshReader reader1 = new MeshReader(gridName+"_omega11.grd");
		
		Mesh meshOmega = reader.read2DMesh();
		Mesh meshOmega0 = reader0.read2DMesh();
		Mesh meshOmega1 = reader1.read2DMesh();
		Mesh meshOmega2 = reader2.read2DMesh();
		
		//Geometry relationship
		meshOmega.computeNodeBelongsToElements();
		meshOmega0.computeNodeBelongsToElements();
		meshOmega1.computeNodeBelongsToElements();
		meshOmega2.computeNodeBelongsToElements();
		meshOmega2.computeNeighborNodes();
		
		//提取外问题meshOmega1的内边界结点的“参数坐标”
		//如果只是提取参数坐标，运行到这里，将输出复制到matlab文件
		//mouse_grid_border.m中的二维数组rs中，运行
		//将matlab输出文件复制到当前目录
		extractBorderParamData(meshOmega1);
		
		//Use element library to assign degree of freedom (DOF) to element
		FELinearTriangleOld triEle = new FELinearTriangleOld();
		//FEBilinearRectangle rectEle = new FEBilinearRectangle();
		FEBilinearRectangleRegular rectEle = new FEBilinearRectangleRegular();
		ElementList eList0 = meshOmega0.getElementList();
		for(int i=1;i<=eList0.size();i++) {
			if(eList0.at(i).nodes.size() == 3)
				triEle.assignTo(eList0.at(i));
			else
				rectEle.assignTo(eList0.at(i));
		}
		ElementList eList2 = meshOmega2.getElementList();
		for(int i=1;i<=eList2.size();i++) {
			if(eList2.at(i).nodes.size() == 3)
				triEle.assignTo(eList2.at(i));
			else
				rectEle.assignTo(eList2.at(i));
		}
		ElementList eList1 = meshOmega1.getElementList();
		for(int i=1;i<=eList1.size();i++) {
			if(eList1.at(i).nodes.size() == 3)
				triEle.assignTo(eList1.at(i));
			else
				rectEle.assignTo(eList1.at(i));
		}
		
		Tools.plotFunction(meshOmega0, outputFolder, "_delta.dat", this.delta);

		//大区域背景解，用于Calibration
		//抽取到Omega上
		Vector u0_bk=null,u_bk=null,u1_bk=null,u2_bk=null;
		//this.mu_a = this.makeMu_a(-0.5, 0.25, 0.3, 0.6, 0.2);
		u0_bk = solveForwardNeumann(meshOmega0);
		Tools.plotVector(meshOmega0, outputFolder, tailType+"_u0_bk.dat", u0_bk);
		//Vector alphaBk0 = this.solveParamInverse(meshOmega0, u0_bk);
		//Tools.plotVector(meshOmega0, outputFolder, "alphaBk0.dat", alphaBk0);
		
		Vector alphaBk =  null;
		if(bUseBaseLineBackground) {
			//-------------求解外问题(Solve exterior problem)-----------
			markExteriorBorder(meshOmega1,tailType);
			//读取测量边界条件，光源编号为srcLightId，时间编号为timeId
			MathFunc diri = null;
			//取timeId=1作为背景
			ObjList<PairDoubleInteger> measBorderValues = readInterpData(srcLightId,1);
			//外问题Omega1的"内外边界"都采用Dirichlet条件
			Vector vDiri = this.extractData(meshOmega0, meshOmega1, u0_bk);
			for(int i=1;i<=measBorderValues.size();i++) {
				PairDoubleInteger di = measBorderValues.at(i);
				vDiri.set(di.i,di.d);
			}
			diri = new Vector2MathFunc(vDiri);
			u1_bk = solveForwardDirichlet(meshOmega1,diri);
			Tools.plotVector(meshOmega1, outputFolder, tailType+"_u1_bk.dat", u1_bk);
			

			u_bk=null;//没有用到，不计算
			
			//注意：由于是求解外问题，Omega2的内部区域（=Omega）没有函数值，
			//不过，后面构造tail的时候，只用到了Omega2的边界值
			u2_bk = this.extractData(meshOmega1, meshOmega2, u1_bk);
			Tools.plotVector(meshOmega2, outputFolder, tailType+"_u2_bk_extract.dat", u2_bk);
		
			//u2_bk通过在Omega2上指定边界条件为测量值，计算获得内部光强值
			HashMap<NodeType, MathFunc> mapNTF = new HashMap<NodeType, MathFunc>();
			mapNTF.put(NodeType.Dirichlet, null);
			meshOmega2.clearBorderNodeMark();
			meshOmega2.markBorderNode(mapNTF);
			diri = new Vector2MathFunc(u2_bk);
			//diri = new Vector2Function(this.extractData(meshOmega0, meshOmega2, u0_bk));//for test
			
			u2_bk = solveForwardDirichlet(meshOmega2,diri);
			Tools.plotVector(meshOmega2, outputFolder, tailType+"_u2_bk.dat", u2_bk);
			//求解系数反问题
			alphaBk = this.solveParamInverse(meshOmega2, u2_bk);
			Tools.plotVector(meshOmega2, outputFolder, "alphaBk.dat", alphaBk);
	
		} else { //从计算出来的u0_bk提取
			//抽取到Omega上
			u_bk = this.extractData(meshOmega0, meshOmega, u0_bk);
			Tools.plotVector(meshOmega, outputFolder, tailType+"_u_bk.dat", u_bk);
			//抽取到Omega1上
			u1_bk = this.extractData(meshOmega0, meshOmega1, u0_bk);
			Tools.plotVector(meshOmega1, outputFolder, tailType+"_u1_bk.dat", u1_bk);
			//抽取到Omega2上
			u2_bk = this.extractData(meshOmega0, meshOmega2, u0_bk);
			Tools.plotVector(meshOmega2, outputFolder, tailType+"_u2_bk.dat", u2_bk);
		}
		
		Vector u1_sim = null;
		if(bSimulate) {
			MathFunc tmp_mu_a = this.mu_a;
			this.mu_a = mu_aSimulate;
			Tools.plotFunction(meshOmega0, outputFolder, tailType+"_mu_a_sim.dat", this.mu_a);
			//大区域模拟Inclusion解，用于观察对比
			Vector u0_sim = solveForwardNeumann(meshOmega0);
			Tools.plotVector(meshOmega0, outputFolder, tailType+"_u0_sim.dat", u0_sim);
			//抽取到Omega上
			Vector u_sim = this.extractData(meshOmega0, meshOmega, u0_sim);
			Tools.plotVector(meshOmega, outputFolder, tailType+"_u_sim.dat", u_sim);
			//抽取到Omega1上
			u1_sim = this.extractData(meshOmega0, meshOmega1, u0_sim);
			Tools.plotVector(meshOmega1, outputFolder, tailType+"_u1_sim.dat", u1_sim);
			//抽取到Omega2上
			Vector u2_sim = this.extractData(meshOmega0, meshOmega2, u0_sim);
			Tools.plotVector(meshOmega2, outputFolder, tailType+"_u2_sim.dat", u2_sim);
			this.mu_a = tmp_mu_a;
		}
		
		
		//------------求解外问题(Solve exterior problem)-------------
		markExteriorBorder(meshOmega1,tailType);
		if(debug) {
			//打印Robin边界条件，检查是否标记正确
			NodeList o1_nodes = meshOmega1.getNodeList();
			for(int i=1;i<=o1_nodes.size();i++) {
				if(o1_nodes.at(i).getNodeType()==NodeType.Robin)
					System.out.println(o1_nodes.at(i));
			}
			//打印Dirichlet边界条件，检查是否标记正确
			for(int i=1;i<=o1_nodes.size();i++) {
				if(o1_nodes.at(i).getNodeType()==NodeType.Dirichlet)
					System.out.println(o1_nodes.at(i));
			}
		}
		//读取测量边界条件，光源编号为srcLightId，时间编号为timeId
		MathFunc diri = null;
		if(bSimulate)
			diri = new Vector2MathFunc(u1_sim);
		else {
			//diri = new DiscreteIndexFunction(readInterpData(srcLightId,timeId));
			ObjList<PairDoubleInteger> measBorderValues = readInterpData(srcLightId,timeId);
			
			//外问题Omega1的"内外边界"都采用Dirichlet条件
			Vector vDiri = u1_bk.copy();
			for(int i=1;i<=measBorderValues.size();i++) {
				PairDoubleInteger di = measBorderValues.at(i);
				vDiri.set(di.i,di.d);
			}
			diri = new Vector2MathFunc(vDiri);
		}
		Vector u1 = solveForwardDirichlet(meshOmega1,diri);
		Tools.plotVector(meshOmega1, outputFolder, tailType+"_u1.dat", u1);
		
		
		//输出Omega1的Dirichlet边界数据
		if(outputMiddleData) {
			NodeList nodesOmega1 = meshOmega1.getNodeList();
			Vector u1Border = new SparseVectorHashMap(u1.getDim());
			for(int i=1;i<=nodesOmega1.size();i++) {
				Node node = nodesOmega1.at(i);
				//if(o1_nodes.at(i).getNodeType()==NodeType.Dirichlet) {
				Variable v = new Variable();
				v.setIndex(i);
				if(tailType == TailType.bottom && isOnBottomBorderOmega(node)) {
					u1Border.set(i,diri.apply(v));
				} else if(tailType == TailType.left && isOnLeftBorderOmega(node)) {
					u1Border.set(i,diri.apply(v));					
				} else if(tailType == TailType.right && isOnRightBorderOmega(node)) {
					u1Border.set(i,diri.apply(v));						
				} else if(tailType == TailType.top && isOnTopBorderOmega(node)) {
					u1Border.set(i,diri.apply(v));						
				}
			}
			Tools.plotVector(meshOmega1, outputFolder, 
					String.format("%s_u1Border%03d.dat", tailType,timeId), u1Border);
		}

		//在规则矩形区域Omega2上，输出用来构造tail的边界上的光强度，测量值延拓后和背景光强
		ObjList<Node> u2BorderNodes = this.getBorderNodes_Omega2(tailType, meshOmega2);
		ObjList<Node> u1BorderNodes = this.getBorderNodes_Omega2(tailType, meshOmega1);
		Vector u2Border = new SparseVectorHashMap(u2_bk.getDim());
		Vector u1Border = new SparseVectorHashMap(u1.getDim());
		for(int i=1;i<=u2BorderNodes.size();i++) {
			int idx = u2BorderNodes.at(i).globalIndex;
			u2Border.set(idx,u2_bk.get(idx));
		}
		for(int i=1;i<=u1BorderNodes.size();i++) {
			int idx = u1BorderNodes.at(i).globalIndex;
			u1Border.set(idx,u1.get(idx));
		}
		Tools.plotVector(meshOmega2, outputFolder, 
				String.format("%s_u2Border_bk%03d.dat", tailType,timeId), u2Border);
		Tools.plotVector(meshOmega1, outputFolder, 
				String.format("%s_u2Border%03d.dat", tailType,timeId), u1Border);
		
		//////////////////////////////////////////////////////////////////////////
	
		//构造tail bottom
		if(tailType==TailType.bottom) {
			Vector tailBottom = computeTailBottom(meshOmega2,u2_bk,meshOmega1,u1);
			//求解系数反问题
			Vector alphaBottom = this.solveParamInverse(meshOmega2, tailBottom);
			if(bUseBaseLineBackground)
				alphaBottom = FMath.axpy(-1.0, alphaBk, alphaBottom);
			Tools.plotVector(meshOmega2, outputFolder, "alphaBottom"+timeId+".dat", alphaBottom);
			
			alphaBottom = Utils.gaussMax(meshOmega2, alphaBottom);
			alphaBottom = Utils.gaussMax(meshOmega2, alphaBottom);
			Tools.plotVector(meshOmega2, outputFolder, "alphaBottom_gaussMax"+timeId+".dat", alphaBottom);
			
			extendAlphaFromBottom(meshOmega2,alphaBottom);
			for(int i=1;i<=alphaBottom.getDim();i++) {
				double v = alphaBottom.get(i);
				if(v<0) alphaBottom.set(i,0);
			}
			Tools.plotVector(meshOmega2, outputFolder, "alphaBottom_gaussMax_extend"+timeId+".dat", alphaBottom);
			tails.add(alphaBottom);
			tailTypes.add("B");
		}
		
		//构造tail top
		if(tailType==TailType.top) {
			Vector tailTop = computeTailTop(meshOmega2,u2_bk,meshOmega1,u1);
			//求解系数反问题
			Vector alphaTop = this.solveParamInverse(meshOmega2, tailTop);
			if(bUseBaseLineBackground)
				alphaTop = FMath.axpy(-1.0, alphaBk, alphaTop);
			Tools.plotVector(meshOmega2, outputFolder, "alphaTop"+timeId+".dat", alphaTop);
			
			alphaTop = Utils.gaussMax(meshOmega2, alphaTop);
			Tools.plotVector(meshOmega2, outputFolder, "alphaTop_gaussMax"+timeId+".dat", alphaTop);
			
			extendAlphaFromBottom(meshOmega2,alphaTop);
			for(int i=1;i<=alphaTop.getDim();i++) {
				double v = alphaTop.get(i);
				if(v<0) alphaTop.set(i,0);
			}
			Tools.plotVector(meshOmega2, outputFolder, "alphaTop_gaussMax_extend"+timeId+".dat", alphaTop);
			tails.add(alphaTop);
			tailTypes.add("T");
		}
		
		//构造tail left
		if(tailType==TailType.left) {
			Vector tailLeft = computeTailLeft(meshOmega2,u2_bk,meshOmega1,u1);
			//求解系数反问题
			Vector alphaLeft = this.solveParamInverse(meshOmega2, tailLeft);
			if(bUseBaseLineBackground)
				alphaLeft = FMath.axpy(-1.0, alphaBk, alphaLeft);
			Tools.plotVector(meshOmega2, outputFolder, "alphaLeft"+timeId+".dat", alphaLeft);
			
			alphaLeft = Utils.gaussMax(meshOmega2, alphaLeft);
			Tools.plotVector(meshOmega2, outputFolder, "alphaLeft_gaussMax"+timeId+".dat", alphaLeft);
			
			extendAlphaFromLeft(meshOmega2,alphaLeft);
			for(int i=1;i<=alphaLeft.getDim();i++) {
				double v = alphaLeft.get(i);
				if(v<0) alphaLeft.set(i,0);
			}
			Tools.plotVector(meshOmega2, outputFolder, "alphaLeft_gaussMax_extend"+timeId+".dat", alphaLeft);
			tails.add(alphaLeft);
			tailTypes.add("L");
		}
		
		//构造tail right
		if(tailType==TailType.right) {
			Vector tailRight = computeTailRight(meshOmega2,u2_bk,meshOmega1,u1);
			//求解系数反问题
			Vector alphaRight = this.solveParamInverse(meshOmega2, tailRight);
			if(bUseBaseLineBackground)
				alphaRight = FMath.axpy(-1.0, alphaBk, alphaRight);
			Tools.plotVector(meshOmega2, outputFolder, "alphaRight"+timeId+".dat", alphaRight);
			
			alphaRight = Utils.gaussMax(meshOmega2, alphaRight);
			Tools.plotVector(meshOmega2, outputFolder, "alphaRight_gaussMax"+timeId+".dat", alphaRight);
			
			extendAlphaFromLeft(meshOmega2,alphaRight);
			for(int i=1;i<=alphaRight.getDim();i++) {
				double v = alphaRight.get(i);
				if(v<0) alphaRight.set(i,0);
			}
			Tools.plotVector(meshOmega2, outputFolder, "alphaRight_gaussMax_extend"+timeId+".dat", alphaRight);
			tails.add(alphaRight);
			tailTypes.add("R");
		}
		if(tails.size() == 2) {
			Vector v1 = tails.at(1);//top,bottom
			Vector v2 = tails.at(2);//left,right
			NodeList nodes2 = meshOmega2.getNodeList();
			//Cut 四周数据
			for(int i=1;i<=v1.getDim();i++) {
				if(nodes2.at(i).coord(1)<-0.75 || nodes2.at(i).coord(1)>0.75) {
					v1.set(i,0);
					v2.set(i,0);
				}
				if(nodes2.at(i).coord(2)<-0.3 || nodes2.at(i).coord(2)>0.7) {
					v1.set(i,0);
					v2.set(i,0);
				}
			}
			//Scale 数据
			double max1 = FMath.max(v1);
			double max2 = FMath.max(v2);
			double min1 = FMath.min(v1);
			double min2 = FMath.min(v2);
			
			double max = Math.max(max1, max2);
			double min = Math.min(min1, min2);
			for(int i=1;i<=v1.getDim();i++) {
				v1.set(i,min+(max-min)*(v1.get(i)-min1)/(max1-min1));
				v2.set(i,min+(max-min)*(v2.get(i)-min2)/(max2-min2));
			}
			
			String sTailTypeTime = String.format("_%s%s%03d", 
					tailTypes.at(1),tailTypes.at(2),timeId);

			if(outputMiddleData) {
				Tools.plotVector(meshOmega2, outputFolder, "alpha"+sTailTypeTime+"_v1.dat", v1);
				Tools.plotVector(meshOmega2, outputFolder, "alpha"+sTailTypeTime+"_v2.dat", v2);
			}
			
			v1 = Utils.gaussSmooth(meshOmega2, v1, 1, 0.5);
			v2 = Utils.gaussSmooth(meshOmega2, v2, 1, 0.5);
			v1.add(1.0, v2);
			
			if(outputMiddleData)
				Tools.plotVector(meshOmega2, outputFolder, "alpha"+sTailTypeTime+"_add.dat", v1);

			v1 = Utils.gaussSmooth(meshOmega2, v1, 1, 0.5);
			
			if(outputMiddleData) {
				Tools.plotVector(meshOmega2, outputFolder, "alpha"+sTailTypeTime+"_add_smooth.dat", v1);
			}
			Tools.plotVector(meshOmega, outputFolder, "final_alpha_omega"+sTailTypeTime+".dat",
					extractData(meshOmega2, meshOmega, v1));
			
			////////////////////////////////////////////////////////////////
//第二种方法：e.g. left, bottom 的外问题解在边界上的值分别延拓的Omega2上后相加，然后求解反问题
//缺点：反问题解不准确			
//			Vector tv1 = tails2.at(1);
//			Vector tv2 = tails2.at(2);
//			for(int i=1;i<=tv1.getDim();i++) {
//				if(tv1.get(i)>0) tv1.set(i,0);
//				if(tv2.get(i)>0) tv2.set(i,0);
//			}
//			max1 = FMath.max(tv1);
//			max2 = FMath.max(tv2);
//			min1 = FMath.min(tv1);
//			min2 = FMath.min(tv2);
//			for(int i=1;i<=tv1.getDim();i++) {
//				tv1.set(i,min2+(max2-min2)*(tv1.get(i)-min1)/(max1-min1));
//				tv2.set(i,min2+(max2-min2)*(tv2.get(i)-min2)/(max2-min2));
//			}
//			Tools.plotVector(meshOmega2, outputFolder, sType+"_tail_tv1.dat", tv1);
//			Tools.plotVector(meshOmega2, outputFolder, sType+"_tail_tv2.dat", tv2);
//			tv1.add(tv2);
//			max1 = FMath.max(tv1);
//			tv1.scale(0.5);
//			tv1.shift(max1);
//			Tools.plotVector(meshOmega2, outputFolder, sType+"_tail_add.dat", tv1);
//			tv1 = Utils.gaussSmooth(meshOmega2, tv1, 2, 0.5);
//			tv1 = Utils.gaussSmooth(meshOmega2, tv1, 2, 0.4);
//			tv1 = Utils.gaussSmooth(meshOmega2, tv1, 2, 0.3);
//			tv1 = Utils.gaussSmooth(meshOmega2, tv1, 2, 0.3);
//			tv1 = Utils.gaussSmooth(meshOmega2, tv1, 2, 0.3);
//			tv1 = Utils.gaussSmooth(meshOmega2, tv1, 2, 0.3);
//			Tools.plotVector(meshOmega2, outputFolder, sType+"_tail_add_smooth.dat", tv1);
//			for(int i=1;i<=nodes2.size();i++) {
//				Node node  = nodes2.at(i);
//				tv1.add(node.globalIndex,u2_bk.get(node.globalIndex));
//			}
//			Tools.plotVector(meshOmega2, outputFolder, sType+"_tail.dat", tv1);
//			//求解系数反问题
//			Vector alpha2 = this.solveParamInverse(meshOmega2, tv1);
//			Tools.plotVector(meshOmega2, outputFolder, "alpha2"+sType+".dat", alpha2);
//			alpha2 = Utils.gaussSmooth(meshOmega2, alpha2, 2, 0.5);
//			alpha2 = Utils.gaussSmooth(meshOmega2, alpha2, 2, 0.4);
//			alpha2 = Utils.gaussSmooth(meshOmega2, alpha2, 2, 0.3);
//			alpha2 = Utils.gaussSmooth(meshOmega2, alpha2, 2, 0.3);
//			alpha2 = Utils.gaussSmooth(meshOmega2, alpha2, 2, 0.3);
//			alpha2 = Utils.gaussSmooth(meshOmega2, alpha2, 2, 0.3);
//			Tools.plotVector(meshOmega2, outputFolder, "alpha2_smooth"+sType+".dat", alpha2);
//			Tools.plotVector(meshOmega, outputFolder, "alpha2_smooth_omega"+sType+".dat",
//					extractData(meshOmega2, meshOmega, alpha2));
			///////////////////////////////////////////////////////////////
		
			tails.clear();
			tails2.clear();
			tailTypes.clear();
		}

	}
	
	public void test() {
		//Read a triangle mesh from an input file
		//包围老鼠头的矩形规则区域
		MeshReader reader2 = new MeshReader("mouse_omega2.grd");
		MeshReader reader00 = new MeshReader("mouse_omega00.grd");
		
		Mesh meshOmega2 = reader2.read2DMesh();
		Mesh meshOmega0 = reader00.read2DMesh();
		
		//Geometry relationship
		meshOmega2.computeNodeBelongsToElements();
		meshOmega2.computeNeighborNodes();
		meshOmega0.computeNodeBelongsToElements();
		meshOmega0.computeNeighborNodes();
		
		//Use element library to assign degree of freedom (DOF) to element
		FELinearTriangleOld linearTriangle = new FELinearTriangleOld();
		ElementList eList2 = meshOmega2.getElementList();
		for(int i=1;i<=eList2.size();i++)
			linearTriangle.assignTo(eList2.at(i));
		ElementList eList0 = meshOmega0.getElementList();
		for(int i=1;i<=eList0.size();i++)
			linearTriangle.assignTo(eList0.at(i));
		
		mu_a=FC.c(0.1);
		Vector u2_bk = solveForwardNeumann(meshOmega2);
		Tools.plotVector(meshOmega2, outputFolder, "test_u2_bk.dat", u2_bk);
		
		mu_a = makeBandMu_a(-0.5,0.1,1,0.1);
		Vector u2 = solveForwardNeumann(meshOmega2);
		Tools.plotVector(meshOmega2, outputFolder, "test_u2.dat", u2);
		Vector u20 = solveForwardNeumann(meshOmega0);
		Tools.plotVector(meshOmega0, outputFolder, "test_u20.dat", u20);
		Vector u2ex = this.extractData(meshOmega0, meshOmega2, u20);
		
		//求解系数反问题
		Vector alphaBottom = this.solveParamInverse(meshOmega2, u2);
		Vector alphaBottomEx = this.solveParamInverse(meshOmega2, u2ex);
		Tools.plotVector(meshOmega2, outputFolder, "test_alphaBottom.dat", alphaBottom);
		Tools.plotVector(meshOmega2, outputFolder, "test_alphaBottomEx.dat", alphaBottomEx);
		
		Tools.plotVector(meshOmega2, outputFolder, "test_diff.dat", 
				u2.add(-1.0, u2_bk));
		
		Tools.plotVector(meshOmega2, outputFolder, "test_diff_ln.dat", 
				FMath.log(FMath.abs(u2)));
	}
	
	/**
	 * 区域说明：
	 * 
	 * Omega: 见PPT文件(D:\RESEARCH\data\2011-3-28Mouse)
	 *        区域范围： x:[-1.5,1.5], y:[-0.4,0.8]
	 * 
	 *     S1--S12--S6
	 *     /          \
	 *    S2           S7
	 *   /              \
	 *  /                \
	 * S4-------S15------r2
	 * 
	 * Omega2:  [-1.6,1.6]*[-0.5,0.9]
	 * Omega00: [-2.1,2.1]*[-1.0,1.4]
	 * Omega11: Omega00\Omega
	 */
	public static void runRat1() {
		MouseHead m = new MouseHead();
		MouseHead.outputMiddleData = true;
		
		//////////////////////Configure/////////////////////////
		//mouse, mouse2 非结构网格 
		//mouse3, mouse4, mouse5 结构网格
		gridName = "mouse4"; 
		
		//760nm波长数据
		m.borderDataPath = "Omega11_Param_T37Groups\\760nm_"+gridName;
		//830nm波长数据
		//m.borderDataPath = "Omega11_Param_T37Groups\\830nm_"+gridName;
		
		m.factor = 10000;
		/////////////////////////////////////////////////////
		
//		m.setDelta(0.0, 0.9);
//		m.test();
		
		//back ground
		m.mu_a = FC.c(0.1);
//		//m.setMu_a(-0.5,0.2,0.1,1);
	
		for(int timeId=1;timeId<=37;timeId++) {
			//bottom tail
//			m.setDelta(-0.9, 0.9); //S1 y(0.8->0.9)
//			m.run(1,timeId,TailType.bottom);
//			m.setDelta(-0.4, 0.9); //S11 y(0.8->0.9)
//			m.run(11,timeId,TailType.bottom);
			
//			m.setDelta(0.0, 0.9);  //S12 y(0.8->0.9)
//			m.run(12,timeId,TailType.bottom);
			
//			m.setDelta(0.4, 0.9);  //S13 y(0.8->0.9)
//			m.run(13,timeId,TailType.bottom);
//			m.setDelta(0.9, 0.9);  //S6 y(0.8->0.9)
//			m.run(6,timeId,TailType.bottom);
			
			//m.setDelta(1.2, 0.4);  //S7 x(1.1->1.2)
			//m.run(7,timeId,TailType.bottom);
			//m.setDelta(-1.2, 0.4); //S2 x(-1.1->-1.2)
			//m.run(2,timeId,TailType.bottom);

			//left tail
//			m.setDelta(1.0, 0.9);  //S6 x(0.9->1.0), y(0.8->0.9)
//			m.run(6,timeId,TailType.left);
			
			m.setDelta(1.2, 0.4);  //S7 x(1.1->1.2)
			m.run(7,timeId,TailType.left);
			
//			m.setDelta(1.4, 0.0);  //S8 x(1.3->1.4)
//			m.run(8,timeId,TailType.left);
//			m.setDelta(1.6, -0.3); //S9 x(1.5->1.6), y(-0.4->-0.3)
//			m.run(9,timeId,TailType.left);

			//right tail
//			m.setDelta(-1.2, 0.4); //S2 x(-1.1->-1.2)
//			m.run(2,timeId,TailType.right);
		}	
	
//		for(int timeId=1;timeId<=4;timeId++) {
//			m.setDelta(0.0, -0.5);//S15 y(-0.4->-0.5)
//			m.run(15,timeId,TailType.top);
//			m.setDelta(-1.2, 0.4);//S2 x(-1.1->-1.2)
//			m.run(2,timeId,TailType.right);
//		}
		
//		m.setDelta(0.8, -0.5);//S16 x(1.3->1.4)
//		m.run(16,1,TailType.left);

	}
	
	public static void runRat6() {
		MouseHead m = new MouseHead();
		MouseHead.outputMiddleData = true;
		
		//////////////////////Configure/////////////////////////
		//mouse, mouse2 非结构网格 
		//mouse3, mouse4, mouse5 结构网格
		gridName = "mouse4"; 
		
		//760nm波长数据 Rat6
		//m.borderDataPath = "Omega11_Param_T31Groups_Rat6\\760nm_"+gridName;
		//830nm波长数据 Rat6
		m.borderDataPath = "Omega11_Param_T31Groups_Rat6\\830nm_"+gridName;

		m.factor = 10000;
		
		/////////////////////////////////////////////////////
		//back ground
		m.mu_a = FC.c(0.1);
		
		for(int timeId=1;timeId<=31;timeId++) {
			//bottom tail
			m.setDelta(0.0, 0.9);  //S13 y(0.8->0.9)
			m.run(13,timeId,TailType.bottom);
//			m.setDelta(-0.4, 0.9);  //S12 y(0.8->0.9)
//			m.run(12,timeId,TailType.bottom);

			//left tail
			m.setDelta(1.2, 0.4);  //S9 x(1.1->1.2)
			m.run(9,timeId,TailType.left);
//			m.setDelta(1.4, 0.0);  //S9 x(1.3->1.4)
//			m.run(10,timeId,TailType.left);
		}
	}
	
	public static void runRat5() {
		MouseHead m = new MouseHead();
		MouseHead.outputMiddleData = true;
		
		//////////////////////Configure/////////////////////////
		//mouse, mouse2 非结构网格 
		//mouse3, mouse4, mouse5 结构网格
		gridName = "mouse4"; 
		
		//760nm波长数据 Rat5
		m.borderDataPath = "Omega11_Param_T21Groups_Rat5\\760nm_"+gridName;
		//830nm波长数据 Rat5
		//m.borderDataPath = "Omega11_Param_T21Groups_Rat5\\830nm_"+gridName;

		m.factor = 10000;
		
		/////////////////////////////////////////////////////
		//back ground
		m.mu_a = FC.c(0.1);
		
		for(int timeId=1;timeId<=21;timeId++) {
			//bottom tail
			m.setDelta(0.0, 0.9);  //S13 y(0.8->0.9)
			m.run(13,timeId,TailType.bottom);
//			m.setDelta(-0.4, 0.9);  //S12 y(0.8->0.9)
//			m.run(12,timeId,TailType.bottom);

			//left tail
			m.setDelta(1.2, 0.4);  //S9 x(1.1->1.2)
			m.run(9,timeId,TailType.left);
//			m.setDelta(1.4, 0.0);  //S9 x(1.3->1.4)
//			m.run(10,timeId,TailType.left);
		}
	}
	public static void main(String[] args) {
		//runRat1();
		//runRat6();
		runRat5();
	}
}
