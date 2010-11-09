package edu.uta.futureye.core;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import edu.uta.futureye.algebra.Vector;
import edu.uta.futureye.core.NodeType;
import edu.uta.futureye.core.intf.GeoEntity;
import edu.uta.futureye.core.intf.Point;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.function.intf.FunctionDerivable;
import edu.uta.futureye.util.Constant;
import edu.uta.futureye.util.DOFList;
import edu.uta.futureye.util.EdgeList;
import edu.uta.futureye.util.ElementList;
import edu.uta.futureye.util.NodeList;
import edu.uta.futureye.util.Utils;

/**
 * 单元
 * @author liuyueming
 *
 */
public class Element {
	/**
	 * Global index of this element
	 */
	public int globalIndex = 0;
	
	/**
	 * Node list of this element
	 * 单元结点列表
	 */
	public NodeList nodes = new NodeList();
	
	/**
	 * Vertices local index in this element(index in nodes)
	 * 单元顶点在nodes中的局部编号列表
	 */
	public List<Integer> verticesLocalIndex = new ArrayList<Integer>(4);
	
	/**
	 * Neighbors of the element
	 * 相邻单元
	 */
	public ElementList neighbors = new ElementList();
	


	/**
	 * A standard method to add a node to the element
	 * 构造一个Element的标准方法：添加一个Node
	 * @param node
	 * @param isVertex
	 */
	public void addNode(Node node, boolean isVertex) {
		nodes.add(node);
		if(isVertex)
			verticesLocalIndex.add(nodes.size());
	}
	
	/**
	 * A map between local index of a node and it's corresponding DOFList
	 * 与单元结点对应的DOFList
	 */
	protected Map<Integer,DOFList> mapDOFList = new LinkedHashMap<Integer,DOFList>();
	
	/**
	 * 来源几何实体
	 */
	public GeoEntity stemGeoEntity = null;
	
	/**
	 * 获取单元结点对应的自由度DOFList
	 * @param localNodeIndex
	 * @return
	 */
	public DOFList getDOFList(int localNodeIndex) {
		DOFList dofList = mapDOFList.get(localNodeIndex);
		return dofList;
	}
	
	/**
	 * 添加一个自由度
	 * @param localNodeIndex
	 * @param dof
	 */
	public void addDOF(int localNodeIndex,DOF dof) {
		DOFList dofList = mapDOFList.get(localNodeIndex);
		if(dofList == null) {
			dofList = new DOFList();
			mapDOFList.put(localNodeIndex, dofList);
		}
		//2010-10-11 DOF反向索引Node
		dof.setOwnerNode(this.nodes.at(localNodeIndex));
		dofList.add(dof);
	}
	
	public void clearDOF() {
		mapDOFList.clear();
	}
	
	
	/**
	 * 局部自由度编号与全局自由度编号转换
	 * @param local
	 * @return
	 */
	public int local2GlobalDOFIndex(int local) {
		for(int j=1;j<=nodes.size();j++) {
			DOFList list = mapDOFList.get(j);
			if(list != null) {
		 		for(int i=1;i<=list.size();i++) {
		 			DOF dof = list.at(i);
					if(dof.getLocalNumber() == local)
						return dof.getGlobalNumber();
				}
	 		}
		}
		return 0;
	}
	
	/**
	 * 获取自由度总数
	 * @return
	 */
	public int getTotalNumberOfDOF() {
		int dim = 0;
		int nNode = nodes.size();
		for(int i=1;i<=nNode;i++) {
			DOFList list = mapDOFList.get(i);
			int nNodeDOF = list.size();
			dim += nNodeDOF;
		}
		return dim;
	}
	
	/**
	 * Return true if exists an edge of this element that is on the border of domain
	 * 判断是否边界单元，即至少存在单元的一边位于区域边界上
	 * @return
	 */
	public boolean isBorderElement() {
		EdgeList edgeList = this.getEdgeList();
		for(int i=1;i<=edgeList.size();i++) {
			Edge edge = edgeList.at(i);
			if(edge.isBorderEdge())
				return true;
		}
		return false;
 	}
	
	public NodeList getNodesByType(NodeType nodeType) {
		NodeList l = new NodeList();
		for(int i=1;i<=nodes.size();i++) {
			if(nodes.at(i).getNodeType() == nodeType) {
				l.add(nodes.at(i));
			}
		}
		return l;
	}
	
	public List<Vertex> getVerticesByType(NodeType nodeType) {
		List<Vertex> l = new LinkedList<Vertex>();
		for(int i=0;i<verticesLocalIndex.size();i++) {
			Node n = nodes.at(verticesLocalIndex.get(i));
			if(n.getNodeType() == nodeType) {
				Vertex v = new Vertex(n.dim());
				v.owner = this;
				v.set(verticesLocalIndex.get(i), n);
				l.add(v);
			}
		}
		return l;
	}
	
	/**
	 * 获取单元的几何顶点，用于坐标变换
	 * @return
	 */
	public List<Vertex> getVertexList() {
		List<Vertex> r = new LinkedList<Vertex>();
		for(Integer index : verticesLocalIndex) {
			Node n = nodes.at(index);
			Vertex v = new Vertex(n.dim());
			//TODO v.owner = this;
			v.set(index, n);
			r.add(v);
		}
		return r;
	}
	
	/**
	 * Get element edges list, dependent on the following
	 * element numbering rules:
	 * 
	 * 4----3
	 * |    |
	 * |    |
	 * 1----2
	 * 
	 * 4--7--3
	 * |     |
	 * 8     6
	 * |     |
	 * 1--5--2
	 * 
	 * 4--11--7--3
	 * |         |
	 * 8         10
	 * |         |
	 * 12        6
	 * |         |
	 * 1-- 5--9--2
	 * 
	 * @return
	 */	
	public EdgeList getEdgeList() {
		EdgeList edgeList = new EdgeList();
		
		List<Vertex> vl = getVertexList();
		//二维单元：
		if(vl.size()>=3 && vl.get(0).dim()==2) {
			vl.add(vl.get(0)); //将第一个顶点再排到最后
			for(int i=1;i<vl.size();i++) {
				Edge edge = new Edge();
				edge.owner = this;
				
				int nodeIndex1 = vl.get(i-1).localIndex;
				int nodeIndex2 = vl.get(i).localIndex;
				
				edge.nodeLocalList.add(
						new NodeLocal(this.nodes.at(nodeIndex1),1)
						);
				edge.nodeLocalList.add(
						new NodeLocal(this.nodes.at(nodeIndex2),2)
				);
				
				int nodeLocalIndex = 3;
				if(vl.size() < nodes.size()) {
					int mul = nodes.size()/vl.size();
					for(int k=0;k<mul;k++) {
						int nodeIndexNN = nodeIndex1 + (k+1)*(vl.size()-1);
						edge.nodeLocalList.add(
							new NodeLocal(this.nodes.at(nodeIndexNN),nodeLocalIndex++)
							);
					}
				}
				edgeList.add(edge);
			}
			return edgeList;
		}
		return null;
	}
	
	public Face getFaceList() {
		return null;
	}
	
	/**
	 * 计算以node为顶点的夹角角度
	 * @param node
	 * @return
	 */
	public double getAngleInElement2D(Node node) {
		int li = getLocalIndex(node);
		int vn = verticesLocalIndex.size();
		if(li <= vn) {
			Node l = nodes.at(li-1<1?vn:li-1);
			Node r = nodes.at(li+1>vn?1:li+1);
			return Utils.computeAngle2D(l, node, r, node);
		} else if(this.nodes.size()/vn == 2){
			Node l = nodes.at(li - vn);
			//TODO 错误的：nodes.at(li - vn + 1)
			Node r = nodes.at( (li - vn + 1)>vn?1:(li - vn + 1));
			return Utils.computeAngle2D(l, node, r, node);
		} else {
			return 0.0; //TODO
		}
		
	}
	
	protected Function jac = null;
	public void updateJacobinLinear1D() {
		String[] fromVars = {"x","y"};
		String[] toVars = {"r"};
		//Coordinate transform and Jacbian on this border element
		CoordinateTransform transBorder = new CoordinateTransform(fromVars,toVars);
		transBorder.transformLinear1D(this);
		jac = (Function) transBorder.getJacobian1D();
	}
	
	public void updateJacobinLinear2D() {
		//Coordinate transform and Jacbian on this element
		CoordinateTransform trans = new CoordinateTransform(2);
		trans.transformLinear2D(this);
		jac = (Function) trans.getJacobian2D();

//TODO adaptive的时候不适用		
//		List<FunctionDerivable> funs = trans.getTransformFunction(
//				trans.getTransformShapeFunctionByElement(this)
//					);
//		trans.setTransformFunction(funs);
//		jac = trans.getJacobian2D();
	}
	
	public Function getJacobin() {
		return jac;
	}
	
	/**
	 * 只处理边界
	 * @param nodeType
	 * @return
	 */
	public ElementList getSubElements() {
		ElementList el = new ElementList();
		EdgeList edgeList = this.getEdgeList();
		for(int i=1;i<=edgeList.size();i++) {
			Edge edge = edgeList.at(i);
			if(edge.isBorderEdge()) {
				Element e = edge.changeToElement();
				el.add(e);
			}
		}
		return el;
	}	
	
	
	public int getLocalIndex(Node node) {
		for(int i=1;i<=nodes.size();i++) {
			if(node.equals(nodes.at(i)))
				return i;
		}
		return 0;
	}
	
	public boolean isVertex(Node node) {
		int localIndex = this.getLocalIndex(node);
		if(localIndex > 0) {
			for(int j=0;j<verticesLocalIndex.size();j++)
				if(verticesLocalIndex.get(j) == localIndex)
					return true;
		}
		return false;
	}
	
	public boolean isBelongElement(Node node) {
		return this.getLocalIndex(node)>0;
	}
	
	public boolean isCoordInElement(Node node) {
		List<Vertex> vList = this.getVertexList();
		//定义循环数组：{1, 2, 3, 1}
		int[] cyc = new int[vList.size()+1];
		for(int k=0;k<cyc.length;k++)
			cyc[k] = k;
		cyc[cyc.length-1] = 0;
		//计算以node为顶点，分别以单元顶点为方向的夹角，如果总和为360，则是内点。
		double angle = 0.0;
		for(int j=0;j<cyc.length-1;j++) {
			Point n1 = vList.get(cyc[j]);
			Point n2 = vList.get(cyc[j+1]);
			angle += Utils.computeAngle2D(n1, node, n2, node);
		}
		if(Math.abs(angle-Math.PI*2)<Constant.eps)
			return true;
		return false;

	}
	public boolean isCoordInElement(double[] coord) {
		Node node = new Node(2);
		node.set(0, coord);
		return isCoordInElement(node);
	}
	
	public String toString() {
		String s = "GE";
		if(globalIndex > 0)
			s += globalIndex;
		s += "( ";
		for(int i=1;i<=nodes.size();i++) {
			String st = "_";
			NodeType nodeType = nodes.at(i).getNodeType();
			if(nodeType == NodeType.Inner)
				st = "I";
			else if(nodeType == NodeType.Dirichlet)
				st = "D";
			else if(nodeType == NodeType.Neumann)
				st = "N";
			else if(nodeType == NodeType.Robin)
				st = "R";
			s += nodes.at(i).globalIndex + st + " ";
		}
		return s+")";
	}
	
	public void adjustVerticeToCounterClockwise() {
		List<Vertex> list = getVertexList();
		if(list.size() <=2) return;
		Vertex v1 = list.get(0);
		Vertex v2 = list.get(1);
		Vertex v3 = list.get(2);
		Vector v12 = null, v13 = null, cp = null;
		if(v1.dim == 2) {
			v12 = Vector.createVector(
					v2.coord(1)-v1.coord(1), v2.coord(2)-v1.coord(2), 0.0);
			v13 = Vector.createVector(
					v3.coord(1)-v1.coord(1), v3.coord(2)-v1.coord(2), 0.0);

		} else if(v1.dim == 3) {
			v12 = Vector.createVector(
					v2.coord(1)-v1.coord(1), v2.coord(2)-v1.coord(2), v2.coord(3)-v1.coord(3)
					);
			v13 = Vector.createVector(
					v3.coord(1)-v1.coord(1), v3.coord(2)-v1.coord(2), v3.coord(3)-v1.coord(3)
					);			
		}
		cp = Vector.crossProduct3D(v12, v13);
		if(cp.get(3) < 0.0) {
			List<Integer> tmp = new ArrayList<Integer>(4);
			//System.out.print(this.toString()+":  ");
			for(int i=verticesLocalIndex.size()-1;i>=0;i--) {
				tmp.add(verticesLocalIndex.get(i));
				//System.out.print(verticesLocalIndex.get(i)+" ");
			}
			//System.out.println("");
			verticesLocalIndex = tmp;
		}
	}
	
	public void addNeighborElement(Element nb) {
		for(int i=1;i<=this.neighbors.size();i++) {
			//TODO ??? nb.globalIndex ???
			if(nb.equals(this.neighbors.at(i)))
				return;
		}
		this.neighbors.add(nb);
	}	
	
	
	
	////////////////////////////////////////////////////////////////////
	/**
	 * 自适应网格加密用来保存从该单元加密出来的子网格单元
	 */
	public ElementList chlids = null;
	public Element parent = null;
	//加密层次
	protected int level = 1;
	
	
	/**
	 * 判断单元是否加密
	 * @return
	 */
	public boolean isRefined() {
		return this.chlids != null;
	}
	
	public int getLevel() {
		return this.level;
	}
	
	public void setLevel(int level) {
		this.level = level;
	}
	////////////////////////////////////////////////////////////////////
	
}
