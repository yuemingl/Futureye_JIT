/**
 * Copyright (c) 2010, nkliuyueming@gmail.com. All rights reserved.
 * 
 * 
 */
package edu.uta.futureye.core;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import edu.uta.futureye.algebra.SparseVectorHashMap;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.tutorial.Tools;
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
	//Mesh name
	protected String name = null;

	//Global node list
	protected NodeList nodeList = new NodeList();
	
	//Global element list
	protected ElementList eleList = new ElementList();
	
	//Global edge list
	protected EdgeList edgeList = null;
	
	//Global face list
	protected FaceList faceList = null;
	
	//Total number of vertices on a mesh
	public int nVertex = 0;
	
	//Boundary types on nodes are defined through functions
	//A function indicates which node belongs to the corresponding NodeType
	//by its coordinates
	protected Map<NodeType, MathFunc> mapNTF;
	
	public boolean debug = false;
	
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
	
	/**
	 * Get a list of all the nodes in the mesh
	 * @return
	 */
	public NodeList getNodeList() {
		return nodeList;
	}
	
	/**
	 * Get a list of all the elements in the mesh
	 * @return
	 */
	public ElementList getElementList() {
		return eleList;
	}
	
	/**
	 * Get a list of boundary nodes in the mesh
	 * @return
	 */
	public NodeList getBoundaryNodeList() {
		NodeList rlt = new NodeList();
		for(int i=1; i<=nodeList.size(); ++i) {
			if(!nodeList.at(i).isInnerNode())
				rlt.add(nodeList.at(i));
		}
		return rlt;
	}
	
	/**
	 * Add a node to the mesh
	 * @param n
	 */
	public void addNode(Node n) {
		nodeList.add(n);
	}	
	/**
	 * Add element <code>e</code> to the mesh, the e.globalIndex will be
	 * assigned automatically
	 *  
	 * @param e
	 */
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
		if(debug)
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
			nodeList.at(i).clearBelongToElements();
		}
		for(int i=1;i<=eleList.size();i++) {
			Element e = eleList.at(i);
			for(int j=1;j<=e.nodes.size();j++) {
				e.nodes.at(j).addBelongToElements(e);
			}
		}
		if(debug)
			System.out.println("computeNodesBelongToElement done!");
	}
	
	public void deleteIsolatedNode() {
		computeNodeBelongsToElements();
		Iterator<Node> node = nodeList.iterator();
		while(node.hasNext()) {
			if(node.next().belongToElements == null) {
				node.remove();
			}
		}
		//Renumbering
		node = nodeList.iterator();
		int index = 1;
		while(node.hasNext()) {
			node.next().globalIndex = index;
			index++;
		}
	}
	
	/**
	 * 添加边界结点类型
	 * @param nodeType
	 * @param fun 控制函数，fun(x)>0边界点，fun(x)<=0内点，x为结点坐标
	 */
	public void addBorderType(NodeType nodeType, MathFunc fun) {
		
	}
	
	public void addBorderType(NodeType nodeType, Node node) {
		
	}
	
	/**
	 * Mark border node type according to mapNTF
	 * 
	 * If the function return a positive value, then the boundary is marked as type indicated in the key of the map
	 * else the boundary is not marked as the type in the key of the map
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
	public void markBorderNode(Map<NodeType,MathFunc> mapNTF) {
		markBorderNode(1,mapNTF);
	}
	
	/**
	 * Mark border node type for component <tt>nVVFComponent</tt> of vector valued unknowns.
	 * 
	 * 对于向量值问题，可以为每个分量<tt>nVVFComponent</tt>分别标记边界类型
	 * 
	 * @param nVVFComponent
	 * @param mapNTF
	 */
	public void markBorderNode(int nVVFComponent, Map<NodeType,MathFunc> mapNTF) {
		if(debug)
			System.out.println("markBorderNode...");
		this.mapNTF = mapNTF;
		if(mapNTF == null) return;
		for(int i=1;i<=nodeList.size();i++) {
			Node node = nodeList.at(i);
			if(node.belongToElements == null) {
				throw new FutureyeException("ERROR: node.belongToElements == null!");
			}
			if(node.belongToElements.size()==0) {
				throw new FutureyeException("ERROR: Call computeNodesBelongToElement() first!");
			} else if(node instanceof NodeRefined && ((NodeRefined) node).isHangingNode()) {
				//TODO Hanging 结点设置为内点，还有其他办法吗？
				node.setNodeType(nVVFComponent, NodeType.Inner); 
			} else {
				if(node.isInnerNode()) {
					node.setNodeType(nVVFComponent, NodeType.Inner); 
				} else { //否则是边界结点
					//System.out.println("Border Node:"+node.globalIndex);
					for(Entry<NodeType,MathFunc> entry : mapNTF.entrySet()) {
						NodeType nodeType = entry.getKey();
						MathFunc fun = entry.getValue();
						if(fun != null) {
							Variable v = new Variable();
							int ic = 1;
							for(String vn : fun.getVarNames())
								v.set(vn, node.coord(ic++));
							if(fun.apply(v) > 0)
								node.setNodeType(nVVFComponent, nodeType);
						} else {
							if(node.getNodeType(nVVFComponent) == null)
								node.setNodeType(nVVFComponent, nodeType);
						}
					}
				}
			}
		}
		if(debug)
			System.out.println("markBorderNode done!");
	}
	
	/**
	 * Mark border node type for components <tt>setVVFComponent</tt> of vector valued unknowns.
	 * 
	 * @param setVVFComponent components indices set
	 * @param mapNTF
	 */
	public void markBorderNode(ObjIndex setVVFComponent, Map<NodeType,MathFunc> mapNTF) {
		for(int i:setVVFComponent)
			markBorderNode(i,mapNTF);
	}
	
	public Map<NodeType, MathFunc> getMarkBorderMap() {
		return this.mapNTF;
	}
	
	public void clearBorderNodeMark() {
		clearBorderNodeMark(1);
	}
	
	public void clearBorderNodeMark(int nVVFComponent) {
		for(int i=1;i<=nodeList.size();i++) {
			Node node = nodeList.at(i);
			if(node.getNodeType() != NodeType.Inner)
				node.setNodeType(nVVFComponent, null);
		}
	}
	public void clearBorderNodeMark(ObjIndex set) {
		for(int i:set)
			clearBorderNodeMark(i);
	}
	
	
	
	/**
	 * Find a node in the mesh that has the same coordination as the given parameter
	 * 
	 * 在网格中按照坐标值查找给定的结点node。如果存在，返回网格中结点对象，否则返回null
	 * 
	 * @param node
	 * @return the node object in the mesh if success, null if no matches
	 */
	public Node findNode(Node node) {
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
	 * Find a node in the mesh that has the same coordination as the given parameter and 
	 * threshold of the distance 
	 * 
	 * 在网格中按照给定的坐标点查找接近的结点。如果存在，返回网格中结点对象，否则返回null
	 * 
	 * @param coord
	 * @param threshold
	 * @return the node object in the mesh if success, null if no matches
	 */
	public Node findNodeByCoord(double[] coord, double threshold) {
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
			nodeList.at(i).clearNeighbors();
		}
		for(int i=1;i<=nodeList.size();i++) {
			Node node = nodeList.at(i);
			if(node.belongToElements==null || node.belongToElements.size()==0) {
				Exception e = new Exception("Call computeNodesBelongToElement() first!");
				e.printStackTrace();
				return;
			}
			else {
				ElementList eList = node.belongToElements;
				for(int j=1;j<=eList.size();j++) {
					Element e = eList.at(j);
					for(int k=1;k<=e.nodes.size();k++) {
						Node nbNode = e.nodes.at(k);
						if(nbNode.globalIndex != node.globalIndex && 
								//2011-11-20 增加是否包含边的判断，可以正确处理四边形单元、六面体单元
								e.containsEdge(node, nbNode))
							node.addNeighbors(nbNode);
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
	 * Copy a mesh
	 */
	public Mesh copy() {
		Mesh newMesh = new Mesh();
		newMesh.nodeList.addAll(nodeList);
		if(edgeList != null) {
			newMesh.edgeList = new EdgeList();
			newMesh.edgeList.addAll(edgeList);
		}
		if(faceList != null) {
			newMesh.faceList = new FaceList();
			newMesh.faceList.addAll(faceList);
		}
		newMesh.eleList.addAll(eleList);
		return newMesh;
	}
	
	/**
	 * Shallow copy or reference
	 * 
	 * @param mesh
	 */
	public void ref(Mesh mesh) {
		this.name = mesh.name;
		this.nodeList = mesh.nodeList;
		this.edgeList = mesh.edgeList;
		this.faceList = mesh.faceList;
		this.eleList = mesh.eleList;
		this.nVertex = mesh.nVertex;
		this.mapNTF = mesh.mapNTF;
		this.debug = mesh.debug;
	}
	
	/**
	 * GEI: Global Element Index
	 * GNI: Global Node Index
	 */
	public void printMeshInfo() {
		for(int i=1;i<=this.eleList.size();i++) {
			Element e = this.eleList.at(i);
			System.out.println("GEI"+i+" "+e);
			for(int j=1;j<=e.nodes.size();j++) {
				Node node = e.nodes.at(j);
				if(node instanceof NodeRefined) {
					NodeRefined nf = (NodeRefined)node;
					if(nf.isHangingNode())
						System.out.println("\t"+nf+" level:"+nf.refineLevel+" HangingNode");
					else
						System.out.println("\t"+nf+" level:"+nf.refineLevel);
						
				} else {
					System.out.println("\t"+node+" level:"+node.refineLevel);
				}
			}
		}
		for(int i=1;i<=this.nodeList.size();i++) {
			System.out.println("GNI"+i+" "+this.nodeList.at(i).globalIndex);
		}
	}
	
	/**
	 * Note:
	 * Inner Node  =5
	 * Dirichlet   =1
	 * Neumann     =2
	 * Robin       =3
	 * Hanging Node=10
	 * Unknown     =20
	 * 
	 * 
	 * @param fileName
	 */
	public void writeNodesInfo(String fileName) {
		writeNodesInfo(fileName,1);
	}
	public void writeNodesInfo(String fileName, int nVVFComponent) {
		int N = nodeList.size();
		Vector v = new SparseVectorHashMap(N,20);
		for(int i=1;i<=N;i++) {
			Node node = nodeList.at(i);
			if(node.getNodeType(nVVFComponent) == NodeType.Inner) 
				v.set(i, 5.0);
			else if(node.getNodeType(nVVFComponent) == NodeType.Dirichlet)
				v.set(i, 1);
			else if(node.getNodeType(nVVFComponent) == NodeType.Neumann)
				v.set(i,2);
			else if(node.getNodeType(nVVFComponent) == NodeType.Robin)
				v.set(i, 3);
			if(node instanceof NodeRefined) {
				NodeRefined nr = (NodeRefined)node;
				if(nr.isHangingNode())
					v.set(i, 10.0);
			}
		}
		Tools.plotVector(this, "", fileName, v);
	}
	
	public String getName() {
		return name;
	}
	
	public void setName(String name) {
		this.name = name;
	}
	
}
