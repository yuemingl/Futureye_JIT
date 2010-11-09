package edu.uta.futureye.core;

import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import edu.uta.futureye.core.intf.Point;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.util.Constant;
import edu.uta.futureye.util.EdgeList;
import edu.uta.futureye.util.ElementList;
import edu.uta.futureye.util.NodeList;
import edu.uta.futureye.util.Utils;

public class Mesh {
	protected NodeList nodeList = new NodeList();
	protected ElementList eleList = new ElementList();
	public static double eps = 1e-6;
	
	public NodeList getNodeList() {
		return nodeList;
	}
	public ElementList getElementList() {
		return eleList;
	}
	public void addNode(Node n) {
		nodeList.add(n);
	}	
	public void addElement(Element e) {
		eleList.add(e);
		e.globalIndex = eleList.size();
	}
	public void clearAll() {
		nodeList.clear();
		eleList.clear();
	}
	
	public void computeNodesBelongToElement() {
		System.out.println("computeNodesBelongToElement...");
		
//		for(int i=1;i<=nodeList.size();i++) {
//			for(int j=1;j<=eleList.size();j++) {
//				Node node = nodeList.at(i);
//				Element e = eleList.at(j);
//				if(e.isInElement(node))
//					node.belongToElements.add(e);
//			}
//			if(i%10==0)
//				System.out.println(String.format("%.1f", 100.0*i/nodeList.size())+"%");
//		}
		
		//New algorithm
		for(int i=1;i<=nodeList.size();i++) {
			nodeList.at(i).belongToElements.clear();
		}
		for(int i=1;i<=eleList.size();i++) {
			Element e = eleList.at(i);
			for(int j=1;j<=e.nodes.size();j++) {
				e.nodes.at(j).addBelongToElements(e);
			}
		}
		
		System.out.println("computeNodesBelongToElement done!");
	}
	
	
	/**
	 * 标记边界结点类型，不会覆盖已标记类型的边界结点，
	 * 如果需要修改，先调用clearBorderNodeMark()
	 * @param mapNTF
	 */
	public void markBorderNode(Map<NodeType,Function> mapNTF) {
		System.out.println("markBorderNode...");
		for(int i=1;i<=nodeList.size();i++) {
			Node node = nodeList.at(i);
			if(node.belongToElements.size()==0)
				System.out.println("ERROR: Call computeNodesBelongToElement() first!");
			else if(node instanceof NodeRefined) {
				//TODO ???
				node.setNodeType(NodeType.Inner); 
			}
			else {
				double sum = 0.0;
				for(int j=1;j<=node.belongToElements.size();j++) {
					sum += node.belongToElements.at(j).getAngleInElement2D(node);
				}
				if(Math.abs(sum-2*Math.PI) <= eps) {
					node.setNodeType(NodeType.Inner); 
				}
				else {
					//System.out.println("Border Node:"+node.globalIndex);
					if(mapNTF != null) {
						for(Entry<NodeType,Function> entry : mapNTF.entrySet()) {
								NodeType nodeType = entry.getKey();
								Function fun = entry.getValue();
							if(fun != null) {
								Variable v = new Variable();
								int ic = 1;
								for(String vn : fun.varNames())
									v.set(vn, node.coord(ic++));
								if(fun.value(v) > 0)
									node.setNodeType(nodeType);
							} else {
								if(node.getNodeType() == null)
									node.setNodeType(nodeType);
							}
						}
					}
				}
			}
		}	
		System.out.println("markBorderNode done!");
	}	
	
	public void clearBorderNodeMark() {
		for(int i=1;i<=nodeList.size();i++) {
			Node node = nodeList.at(i);
			if(node.getNodeType() != NodeType.Inner)
				node.setNodeType(null);
		}
	}
	
	
	/**
	 * 判断网格是否包含结点node
	 * @param node
	 * @return null if not contains the node
	 */
	public Node containNode(Node node) {
		for(int i=1;i<=nodeList.size();i++) {
			int nDim = node.dim();
			boolean same = true;
			for(int j=1;j<=nDim;j++) {
				if(Math.abs(node.coord(j)-nodeList.at(i).coord(j)) > Constant.eps) {
					same = false;
					break;
				}
			}
			if(same) 
				return nodeList.at(i);
		}
		return null;
	}
	
	/**
	 * 获取与该坐标点接近的结点
	 * @param coord
	 * @param threshold
	 * @return
	 */
	public Node getNodeByCoord(double[] coord, double threshold) {
		for(int i=1;i<=nodeList.size();i++) {
			int nDim = nodeList.at(1).dim();
			boolean same = true;
			for(int j=1;j<=nDim;j++) {
				if(Math.abs(coord[j-1]-nodeList.at(i).coord(j)) > threshold) {
					same = false;
					break;
				}
			}
			if(same)
				return nodeList.at(i);
		}
		return null;
	}
	
	/**
	 * 获取包含该坐标点的单元（二维）
	 * @param coord
	 * @return
	 */
	public Element getElementByCoord(double[] coord) {
		Node node = new Node(2);
		node.set(0, coord);
		for(int i=1;i<=eleList.size();i++) {
			Element e = eleList.at(i);
			if(e.isCoordInElement(coord))
				return e;
		}
		return null;
	}
	
	public void computeNeiborNode() {
		for(int i=1;i<=nodeList.size();i++) {
			nodeList.at(i).neighbors.clear();
		}
		for(int i=1;i<=nodeList.size();i++) {
			Node node = nodeList.at(i);
			if(node.belongToElements.size()==0) {
				Exception e = new Exception("Call computeNodesBelongToElement() first!");
				e.printStackTrace();
				return;
			}
			else {
				ElementList eList = node.belongToElements;
				for(int j=1;j<=eList.size();j++) {
					NodeList nList = eList.at(j).nodes;
					for(int k=1;k<=nList.size();k++) {
						Node nbNode = nList.at(k);
						if(nbNode.globalIndex != node.globalIndex)
							node.neighbors.add(nbNode);
					}
				}
			}
		}
	}
	
	public void computeNeighborElement() {
		for(int i=1;i<=eleList.size();i++) {
			eleList.at(i).neighbors.clear();
		}
		for(int i=1;i<=nodeList.size();i++) {
			Node node = nodeList.at(i);
			if(node.belongToElements.size()==0) {
				Exception e = new Exception("Call computeNodesBelongToElement() first!");
				e.printStackTrace();
				return;
			}
			else {
				ElementList eList = node.belongToElements;
				for(int j=1;j<=eList.size();j++) {
					for(int k=1;k<=eList.size();k++) {
						Element e1 = eList.at(j);
						Element e2 = eList.at(k);
						if(e1.equals(e2))
							continue;
						if(isNeighbor(e1,e2)) {
							e1.addNeighborElement(e2);
							e2.addNeighborElement(e1);
						}
					}
				}
			}
		}
		
	}
	
	/**
	 * 仅适用于二维
	 * @param e1
	 * @param e2
	 * @return
	 */
	public boolean isNeighbor(Element e1, Element e2) {
		int counter = 0;
//		for(int i=1;i<=e1.nodes.size();i++) {
//			for(int j=1;j<=e2.nodes.size();j++) {
//				//有两个公共结点(对于自适应网格，由于存在hanging node，它不属于粗单元格，因此该方法判断有误)
//				if(e1.nodes.at(i).globalIndex == e2.nodes.at(j).globalIndex) {
//					counter++;
//					if(counter==2) return true; 
//				}
//			}
//		}
		EdgeList e1EdgeList = e1.getEdgeList();
		EdgeList e2EdgeList = e2.getEdgeList();
		for(int i=1;i<=e1EdgeList.size();i++) {
			for(int j=1;j<=e2EdgeList.size();j++) {
				NodeList e1EndNodes = e1EdgeList.at(i).getEndNodes();
				NodeList e2EndNodes = e2EdgeList.at(j).getEndNodes();
				if(Utils.isLineOverlap(e1EndNodes.at(1), e1EndNodes.at(2), 
						e2EndNodes.at(1), e2EndNodes.at(2)))
					return true;
			}
		}
		return false;
	}
}
