package edu.uta.futureye.application;

import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;

import edu.uta.futureye.algebra.SparseVectorHashMap;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.core.EdgeLocal;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.core.NodeType;
import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.util.Constant;
import edu.uta.futureye.util.FutureyeException;
import edu.uta.futureye.util.PairDoubleInteger;
import edu.uta.futureye.util.Utils;
import edu.uta.futureye.util.container.ElementList;
import edu.uta.futureye.util.container.NodeList;
import edu.uta.futureye.util.container.ObjList;

public class MouseHead20110908 extends MouseHead {
	
	/**
	 * 提取网格Omega1（外问题）内边界点的参数坐标，提供给matlab程序的
	 * “插值结点的参数坐标(0->r1->r2->r3->r4(0))”
	 * 
	 *  S25            S1
	 * (r4)0---->----r1
	 *     /          \
	 *    /            \
	 *   ┐              ┘
	 *  /                \
	 * r3-------<--------r2
	 * S28               S4
	 * 
	 * 在标准输出打印结点编号、参数坐标对列表，以参数坐标值顺序排序
	 */	
	@Override
	public void extractBorderParamData(Mesh omega) {
		ObjList<PairDoubleInteger> rs = new ObjList<PairDoubleInteger>();
		NodeList nodes = omega.getNodeList();
		Node S25 = new Node(0,-1.1, 1.0);
		Node S1  = new Node(0, 1.1, 1.0);
		Node S4  = new Node(0, 1.5,-0.5);
		Node S28 = new Node(0,-1.5,-0.5);
		double dis1 = Utils.computeLength(S25,S1);
		double dis2 = dis1 + Utils.computeLength(S1,S4);
		double dis3 = dis2 + Utils.computeLength(S4,S28);
		//double dis4 = dis3 + Utils.computeLength(S4,S1);
		
		for(int i=1;i<=nodes.size();i++) {
			Node node = nodes.at(i);
			if(node.isInnerNode())
				continue;
			//除去小矩形（外问题在Omega1上求解时）
			if(Math.abs(node.coord(1)-(-1.6))<Constant.meshEps ||
					Math.abs(node.coord(1)-1.6)<Constant.meshEps ||
					Math.abs(node.coord(2)-(-0.6))<Constant.meshEps ||
					Math.abs(node.coord(2)-1.1)<Constant.meshEps)
				continue;
			//除去大矩形（外问题在Omega11上求解）
			if(Math.abs(node.coord(1)-(-2.6))<Constant.meshEps ||
					Math.abs(node.coord(1)-2.6)<Constant.meshEps ||
					Math.abs(node.coord(2)-(-1.6))<Constant.meshEps ||
					Math.abs(node.coord(2)-2.1)<Constant.meshEps)
				continue;
			
			PairDoubleInteger di = new PairDoubleInteger();
			di.i = node.getIndex();
			
			if(Utils.isPointOnLineSegment(S25, S1, node)) {//top
				di.d = Utils.computeLength(S25,node);
			} else if (Utils.isPointOnLineSegment(S1, S4, node)) {//right
				di.d = dis1 + Utils.computeLength(S1,node);
			} else if(Utils.isPointOnLineSegment(S4, S28, node)) {//bottom
				di.d = dis2 + Utils.computeLength(S4,node);
			} else if(Utils.isPointOnLineSegment(S28, S25, node)) {//left
				di.d = dis3 + Utils.computeLength(S28,node);
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
	 * 根据不同的tail标记不同的边界类型
	 * @param mesh
	 * @param fTailType
	 */
	@Override
	public void markExteriorBorder(Mesh mesh, final TailType fTailType) {
		//外问题 Solve exterior problem
		HashMap<NodeType, MathFunc> mapNTF = new HashMap<NodeType, MathFunc>();
//		mapNTF.put(NodeType.Robin, new AbstractFunction("x","y") {
//			@Override
//			public double value(Variable v) {
//				double x = v.get("x");
//				double y = v.get("y");
//				//Omega11 外问题区域的外边界
//				if(Math.abs(x-2.6)<Constant.meshEps ||
//						Math.abs(x+2.6)<Constant.meshEps ||
//						Math.abs(y-2.1)<Constant.meshEps ||
//						Math.abs(y+1.6)<Constant.meshEps)
//					return 1;
//				//bottom tail: 外问题上边界不要Dirichlet条件，改为Robin
//				if(fTailType==TailType.bottom && Math.abs(y-1.0)<Constant.meshEps)
//					return 1;
//				//top tail: 外问题右边界不要Dirichlet条件，改为Robin
//				if(fTailType==TailType.top && Math.abs(y-(-0.5))<Constant.meshEps)
//					return 1;
//				//最后处理left tail和right tail，两条边界是斜线，只能指定一个范围：
//				//left tail: 外问题右边界不要Dirichlet条件，改为Robin
//				if(fTailType==TailType.left && x > 1.1-Constant.meshEps)
//					return 1;
//				//right tail: 外问题右边界不要Dirichlet条件，改为Robin
//				if(fTailType==TailType.right && x < -1.1+Constant.meshEps)
//					return 1;
//				
//				return 0;
//			}
//		});
		
		mapNTF.put(NodeType.Dirichlet, null);
		mesh.clearBorderNodeMark();
		mesh.markBorderNode(mapNTF);		
	}
	
	
	/**
	 * 获取区域的边界结点列表（与Omega2边界一致的任何区域）
	 * Omega2:  [-1.6,1.6]*[-0.5,0.9]
	 * 
	 * @param type
	 * @param omega2
	 * @return
	 */
	@Override
	public ObjList<Node> getBorderNodes_Omega2(TailType type, Mesh omega2) {
		NodeList nodes2 = omega2.getNodeList(); 
		ObjList<Node> borderNodes = new ObjList<Node>();
		double x,y;
		if(type == TailType.left) {
			for(int i=1;i<=nodes2.size();i++) {
				if( //y=[-0.6,1.1]
					nodes2.at(i).coord(2)-(-0.6)>-Constant.meshEps && 1.1-nodes2.at(i).coord(2)>-Constant.meshEps &&
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
					//y:-0.6
					Math.abs(y-(-0.6))<Constant.meshEps
				) {
					borderNodes.add(nodes2.at(i));
					if(debug)
						System.out.println(nodes2.at(i));
				}
			}
		}  else if(type == TailType.right) {
			for(int i=1;i<=nodes2.size();i++) {
				if( //y=[-0.6,1.1]
					nodes2.at(i).coord(2)-(-0.6)>-Constant.meshEps && 1.1-nodes2.at(i).coord(2)>-Constant.meshEps &&
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
					//y:1.1
					Math.abs(nodes2.at(i).coord(2)-(1.1))<Constant.meshEps
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
			
			if(Math.abs(node.coord(2)-(-0.6))<Constant.meshEps) 
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
							
							double MaxS = Math.sqrt((ly-(-0.6))*(ly-(-0.6)));
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
			
			if(Math.abs(node.coord(2)-1.1)<Constant.meshEps) 
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
							
							double MaxS = Math.sqrt((ly-1.1)*(ly-1.1));
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
			
			if(Math.abs(node.coord(2)-(-0.6))<Constant.meshEps) 
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
			
			if(Math.abs(node.coord(2)-1.1)<Constant.meshEps) 
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
	
	public static void runRat1(String waveLength, int nTotalTime) {
		MouseHead20110908 m = new MouseHead20110908();
		MouseHead20110908.outputMiddleData = true;
		MouseHead20110908.debug = true;
		
		//////////////////////Configure/////////////////////////
		//结构网格
		gridName = "mouse0908_1"; 
		
		//waveLength nm 波长数据
		m.borderDataPath = "Omega11_Param_T"+nTotalTime+"Groups_Rat1/"+waveLength+"nm_"+gridName;

		m.outputFolder = "MouseHead20110908";
		m.factors = new double[29];
		
		m.baselines = new double[29];
		m.factors[9] = 20000;
		m.baselines[9] = 0;
		m.factors[17] = 20000;
		m.baselines[17] = 0;
		m.factors[13] = 20000;
		m.baselines[13] = 0;
		
		m.factors[2] = 20000;
		m.baselines[2] = 0;
		
		//back ground of mu_a and mu_s
//for bUseBaseLineBackground = false;		
//		m.mu_a = FC.c(0.1);
//		m.k = new FC(1.0/100.0);
//for bUseBaseLineBackground = true;		
		m.mu_a = FC.c(0.1);
		m.mu_a_bk = 0.1;
		m.k = new FC(1.0/50.0);
		
		m.bSimulate = false; //是否使用模拟数据
		m.mu_aSimulate = m.makeMu_a(-0.5, 0.25, 0.3, 0.6, 0.1);
		//m.mu_aSimulate = m.makeMu_a(-1.0, 0.25, 0.3, 0.6, 0.1);
		
		m.bUseBaseLineBackground = true;
		
		/////////////////////////////////////////////////////
		
		for(int timeId=1;timeId<=nTotalTime;timeId++) {
			//bottom tail
			m.setDelta(0.0, 0.9);  //S13 y(1.1->1.5)
			m.run(13,timeId,TailType.bottom);
//			m.setDelta(-0.5, 0.9);  //S17 y(1.1->1.2)
//			m.run(17,timeId,TailType.bottom);
//			m.setDelta(0.5, 0.9);  //S9 y(1.1->1.2)
//			m.run(9,timeId,TailType.bottom);

			//left tail
			m.setDelta(1.033, 0.5);  //S2 x(1.233->1.3)
			m.run(2,timeId,TailType.left);
//			m.setDelta(1.467, 0.0);  //S2 x(1.367->1.467)
//			m.run(3,timeId,TailType.left);
		}
	}
	
	public static void runRat2(String waveLength, int nTotalTime) {
		MouseHead20110908 m = new MouseHead20110908();
		MouseHead20110908.outputMiddleData = true;
		MouseHead20110908.debug = true;
		
		//////////////////////Configure/////////////////////////
		//结构网格
		gridName = "mouse0908_1"; 
		
		//760nm波长数据 Rat5
		m.borderDataPath = "Omega11_Param_T"+nTotalTime+"Groups_Rat2/"+waveLength+"nm_"+gridName;

		m.outputFolder = "MouseHead20110908";
		m.factors = new double[29];
		
		m.baselines = new double[29];
		m.factors[9] = 20000;
		m.baselines[9] = 0;
		m.factors[17] = 20000;
		m.baselines[17] = 0;
		m.factors[13] = 20000;
		m.baselines[13] = 0;
		
		m.factors[2] = 20000;
		m.baselines[2] = 0;
		
		//back ground of mu_a and mu_s
		//for bUseBaseLineBackground = false;		
//		m.mu_a = FC.c(0.1);
//		m.k = new FC(1.0/100.0);
//for bUseBaseLineBackground = true;
		
//max(mu_a)=0.6		
//		m.mu_a = FC.c(0.22);
//		m.mu_a_bk = 0.22;
//		m.k = new FC(1.0/40.0);
//max(mu_a)=0.3		
		m.mu_a = FC.c(0.1);
		m.mu_a_bk = 0.1;
		m.k = new FC(1.0/50.0);
		
		m.bSimulate = false; //是否使用模拟数据
		m.mu_aSimulate = m.makeMu_a(-0.5, 0.25, 0.3, 0.6, 0.1);
		//m.mu_aSimulate = m.makeMu_a(-1.0, 0.25, 0.3, 0.6, 0.1);
		
		m.bUseBaseLineBackground = true;

		/////////////////////////////////////////////////////
		
		for(int timeId=1;timeId<=nTotalTime;timeId++) {
			//bottom tail
			m.setDelta(0.0, 0.9);  //S13 y(1.1->1.5)
			m.run(13,timeId,TailType.bottom);
//			m.setDelta(-0.5, 0.9);  //S17 y(1.1->1.2)
//			m.run(17,timeId,TailType.bottom);
//			m.setDelta(0.5, 0.9);  //S9 y(1.1->1.2)
//			m.run(9,timeId,TailType.bottom);

			//left tail
			m.setDelta(1.033, 0.5);  //S2 x(1.233->1.3)
			m.run(2,timeId,TailType.left);
//			m.setDelta(1.467, 0.0);  //S2 x(1.367->1.467)
//			m.run(3,timeId,TailType.left);
		}
	}
	
	public boolean isOnTopBorderOmega(Node node) {
		double x = node.coord(1);
		double y = node.coord(2);
		if(Math.abs(y-1.0)<Constant.meshEps &&
				x>=-1.1-Constant.meshEps && x<=1.1+Constant.meshEps)
			return true;
		else
			return false;
	}
	
	public boolean isOnBottomBorderOmega(Node node) {
		double x = node.coord(1);
		double y = node.coord(2);
		if(Math.abs(y-(-0.5))<Constant.meshEps &&
				x>=-1.5-Constant.meshEps && x<=1.5+Constant.meshEps)
			return true;
		else
			return false;
	}
	
	public boolean isOnLeftBorderOmega(Node node) {
		Node S25 = new Node(0,-1.1,1.0);
		Node S28 = new Node(0,-1.5,-0.5);
		return Utils.isPointOnLineSegment(S25, S28, node);
	}
	
	public boolean isOnRightBorderOmega(Node node) {
		Node S1 = new Node(0,1.1,1.0);
		Node S4 = new Node(0,1.5,-0.5);
		return Utils.isPointOnLineSegment(S1, S4, node);
	}	
	
	public static void main(String[] args) {
		if(args.length==0)
			System.out.println("java -jar runMouseHead20110908.jar [rat1|rat2] [760|830]");// [mT-nT]
		else {
			if(args[0].equals("rat1"))
				runRat1(args[1],261);
			else if(args[0].equals("rat2"))
				runRat2(args[1],341);
			else
				System.out.println("args[0]=rat1 or rat2");
		}
	}
}
