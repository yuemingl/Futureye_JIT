package edu.uta.futureye.core;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.Map.Entry;

import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.util.Constant;
import edu.uta.futureye.util.FutureyeException;
import edu.uta.futureye.util.MultiKey;
import edu.uta.futureye.util.Utils;
import edu.uta.futureye.util.container.EdgeList;
import edu.uta.futureye.util.container.ElementList;
import edu.uta.futureye.util.container.FaceList;
import edu.uta.futureye.util.container.NodeList;
import edu.uta.futureye.util.container.ObjIndex;
import edu.uta.futureye.util.container.ObjList;

public class Mesh {
	//Global node list
	protected NodeList nodeList = new NodeList();
	
	//Global element list
	protected ElementList eleList = new ElementList();
	
	//Global edge list
	protected EdgeList edgeList = null;
	
	//Global face list
	protected FaceList faceList = null;
	
	public int nVertex = 0;
	
	public EdgeList getEdgeList() {
		return edgeList;
	}
	public void setEdgeList(EdgeList edgeList) {
		this.edgeList = edgeList;
	}
	
	public FaceList getFacelist() {
		return faceList;
	}
	public void setFacelist(FaceList facelist) {
		this.faceList = facelist;
	}

	//Global volume list
	//??? protected VolumeList volumeList = null;
	
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
	
	/**
	 * 计算结点所属的单元，在计算其他网格关系时，该步骤必须先计算
	 */
	public void computeNodeBelongsToElements() {
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
	 * 添加边界结点类型
	 * @param nodeType
	 * @param fun 控制函数，fun(x)>0边界点，fun(x)<=0内点，x为结点坐标
	 */
	public void addBorderType(NodeType nodeType, Function fun) {
		
	}
	
	public void addBorderType(NodeType nodeType, Node node) {
		
	}
	
	/**
	 * Mark border node type according to mapNTF
	 * 
	 * 标记边界结点类型，不会覆盖已标记过类型的边界结点，
	 * 如果需要修改，先调用clearBorderNodeMark()
	 * 
	 * 依赖：computeNodesBelongToElement()
	 * 
	 * TODO
	 * 如果网格本身包含边界结点类型，如何标记？
	 * 
	 * @param mapNTF
	 */
	public void markBorderNode(Map<NodeType,Function> mapNTF) {
		markBorderNode(1,mapNTF);
	}
	
	/**
	 * 对于向量值问题，可以为每个分量<tt>vvfIndex</tt>分别标记边界类型
	 * 
	 * @param vvfIndex
	 * @param mapNTF
	 */
	public void markBorderNode(int vvfIndex, Map<NodeType,Function> mapNTF) {
		System.out.println("markBorderNode...");
		for(int i=1;i<=nodeList.size();i++) {
			Node node = nodeList.at(i);
			if(node.belongToElements.size()==0) {
				throw new FutureyeException("ERROR: Call computeNodesBelongToElement() first!");
			} else if(node instanceof NodeRefined && ((NodeRefined) node).isHangingNode()) {
				//TODO Hanging 结点设置为内点，还有其他办法吗？
				node.setNodeType(vvfIndex, NodeType.Inner); 
			} else {
				if(node.isInnerNode()) {
					node.setNodeType(vvfIndex, NodeType.Inner); 
				} else { //否则是边界结点
					//System.out.println("Border Node:"+node.globalIndex);
					for(Entry<NodeType,Function> entry : mapNTF.entrySet()) {
						NodeType nodeType = entry.getKey();
						Function fun = entry.getValue();
						if(fun != null) {
							Variable v = new Variable();
							int ic = 1;
							for(String vn : fun.varNames())
								v.set(vn, node.coord(ic++));
							if(fun.value(v) > 0)
								node.setNodeType(vvfIndex, nodeType);
						} else {
							if(node.getNodeType(vvfIndex) == null)
								node.setNodeType(vvfIndex, nodeType);
						}
					}
				}
			}
		}
		System.out.println("markBorderNode done!");
	}	
	public void markBorderNode(ObjIndex vvfIndexSet, Map<NodeType,Function> mapNTF) {
		for(int i:vvfIndexSet)
			markBorderNode(i,mapNTF);
	}
	
	public void clearBorderNodeMark() {
		clearBorderNodeMark(1);
	}
	
	public void clearBorderNodeMark(int vvfIndex) {
		for(int i=1;i<=nodeList.size();i++) {
			Node node = nodeList.at(i);
			if(node.getNodeType() != NodeType.Inner)
				node.setNodeType(vvfIndex, null);
		}
	}
	public void clearBorderNodeMark(ObjIndex set) {
		for(int i:set)
			clearBorderNodeMark(i);
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
				if(Math.abs(node.coord(j)-nodeList.at(i).coord(j)) > Constant.meshEps) {
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
	
	public Element getElementByNodes(NodeList nodes) {
		for(int i=1;i<=this.eleList.size();i++) {
			Element e = this.eleList.at(i);
			boolean find = true;
			if(e.nodes.size()!=nodes.size())
				continue;
			for(int j=1;j<=nodes.size();j++) {
				if(!e.nodes.at(j).coordEquals(nodes.at(j))) {
					find = false;
					break;
				}
			}
			if(find)
				return e;
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
	
	/**
	 * Compute neighbor nodes of a node
	 * 计算相邻结点
	 * 
	 * Depends:
	 *   computeNodesBelongToElement()
	 * 
	 */
	public void computeNeighborNodes() {
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
	
	/**
	 * Compute neighbor elements of an element
	 * 计算相邻单元
	 * 
	 * Depends:
	 *  computeNodesBelongToElement()
	 *  computeGlobalEdge() (2D,3D case)
	 *  computeGlobalFace() (3D case)
	 *
	 */
	public void computeNeighborElements() {
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
	 * Compute global edges in grid
	 * 计算网格包含的全局边
	 * 
	 */
	public void computeGlobalEdge() {
		Map<MultiKey, Edge> map = new HashMap<MultiKey, Edge>();
		for(int i=1;i<=eleList.size();i++) {
			Element e = eleList.at(i);
			ObjList<EdgeLocal> localEdges = e.edges();//单元e的局部边列表
			for(int j=1;j<=localEdges.size();j++) {
				EdgeLocal localEdge = localEdges.at(j);//局部边
				MultiKey mkey = new MultiKey(
								localEdge.getVertices().at(1).globalNode().globalIndex,
								localEdge.getVertices().at(2).globalNode().globalIndex
						);
				Edge globalEdge = map.get(mkey);
				if(globalEdge == null) {
					globalEdge = localEdge.buildEdge();
					map.put(mkey, globalEdge);
				}
				//update global edge
				localEdge.globalEdge = globalEdge;
			}
		}
		//为mesh.edgeList（全局边）赋值，并加全局边索引（编号）
		this.edgeList = new EdgeList();
		int globalIndex = 1;
		for(Entry<MultiKey, Edge> entry : map.entrySet()) {
			//globalIndex：全局边索引（编号）
			//globalIndex: Global index of Global edges
			entry.getValue().setGlobalIndex(globalIndex++);
			this.edgeList.add(entry.getValue());
		}
	}
	
	/**
	 * Compute global faces in grid
	 * 计算网格包含的全局面
	 * 
	 */
	public void computeGlobalFace() {
		Map<Set<Integer>, Face> map = new HashMap<Set<Integer>, Face>();
		for(int i=1;i<=eleList.size();i++) {
			Element e = eleList.at(i);
			//遍历局部面，利用HashMap做face的剔重判断，得到全局面
			ObjList<FaceLocal> localFaces = e.faces();
			for(int j=1; j<=localFaces.size(); j++) {
				FaceLocal localFace = localFaces.at(j);
				Set<Integer> nodeSet = new HashSet<Integer>();
				ObjList<Vertex> vertices = localFace.getVertices();
				for(int n=1;n<=vertices.size();n++) {
					nodeSet.add(vertices.at(n).globalNode().globalIndex);
				}
				Face gobalFace = map.get(nodeSet);
				if(gobalFace == null) {
					gobalFace = localFace.buildFace();
					map.put(nodeSet, gobalFace);
				}
				//update global edge
				localFace.globalFace = gobalFace;
			}
		}
		//为mesh.faceList（全局面）赋值，并加全局面索引（编号）
		this.faceList = new FaceList();
		int globalIndex = 1;
		for(Entry<Set<Integer>, Face> entry : map.entrySet()) {
			//globalIndex：全局面索引（编号）
			//globalIndex: Global index of Global faces
			entry.getValue().setGlobalIndex(globalIndex++);
			this.faceList.add(entry.getValue());
		}
	}
	
	/**
	 * 仅适用于二维
	 * @param e1
	 * @param e2
	 * @return
	 */
	//2011-02-19
	public boolean isNeighbor(Element e1, Element e2) {
//		int counter = 0;
//		for(int i=1;i<=e1.nodes.size();i++) {
//			for(int j=1;j<=e2.nodes.size();j++) {
//				//有两个公共结点(对于自适应网格，由于存在hanging node，它不属于粗单元格，因此该方法判断有误)
//				if(e1.nodes.at(i).globalIndex == e2.nodes.at(j).globalIndex) {
//					counter++;
//					if(counter==2) return true; 
//				}
//			}
//		}
		ObjList<EdgeLocal> e1EdgeList = e1.edges();
		ObjList<EdgeLocal> e2EdgeList = e2.edges();
		for(int i=1;i<=e1EdgeList.size();i++) {
			for(int j=1;j<=e2EdgeList.size();j++) {
				NodeList e1EndNodes = e1EdgeList.at(i).globalEdge.getEndNodes();
				NodeList e2EndNodes = e2EdgeList.at(j).globalEdge.getEndNodes();
				if(Utils.isLineOverlap(e1EndNodes.at(1), e1EndNodes.at(2), 
						e2EndNodes.at(1), e2EndNodes.at(2)))
					return true;
			}
		}
		return false;
	}
	
	/**
	 * 
	 */
	Mesh copy() {
		//TODO
		return null;
	}
	
	public void printMeshInfo() {
		for(int i=1;i<=this.eleList.size();i++) {
			Element e = this.eleList.at(i);
			System.out.println("GEI"+i+" "+e);
			for(int j=1;j<=e.nodes.size();j++) {
				Node node = e.nodes.at(j);
				if(node instanceof NodeRefined) {
					NodeRefined nf = (NodeRefined)node;
					if(nf.isHangingNode())
						System.out.println("\t"+nf+" level:"+nf.level+" HangingNode");
					else
						System.out.println("\t"+nf+" level:"+nf.level);
						
				} else {
					System.out.println("\t"+node+" level:"+node.level);
				}
			}
		}
		for(int i=1;i<=this.nodeList.size();i++) {
			System.out.println("GNI"+i+" "+this.nodeList.at(i).globalIndex);
		}
	}
	
}
