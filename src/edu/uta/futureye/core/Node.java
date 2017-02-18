/**
 * Copyright (c) 2010, nkliuyueming@gmail.com. All rights reserved.
 * 
 * 
 */
package edu.uta.futureye.core;

import edu.uta.futureye.core.geometry.Point;
import edu.uta.futureye.util.Constant;
import edu.uta.futureye.util.FutureyeException;
import edu.uta.futureye.util.container.ElementList;
import edu.uta.futureye.util.container.NodeList;
import edu.uta.futureye.util.container.ObjVector;

/**
 * Global finite element node class
 * 有限元（全局）结点
 * 
 * @author liuyueming
 *
 */
public class Node implements Point {
	//全局索引（编号）
	public int globalIndex = 0;
	
	//全局来看，结点所属单元
	public ElementList belongToElements = null;
	
	//相邻结点
	public NodeList neighbors = null;
	
	//空间维度
	protected int dim = 0;
	
	//空间坐标
	protected double[] coords = new double[3];
	
	//结点类型：边界结点 or 内部结点
	protected ObjVector<NodeType> nodeTypes = null;
	
	//结点属于的加密层次，即第level次加密产生的结点
	protected int refineLevel = 1;

	public Node() {
	}
	
	public Node(int dim) {
		this.dim = dim;
	}
	
	/**
	 * 构造一个全局结点
	 * @param globalIndex 全局索引（编号）
	 * @param x 第一个坐标
	 * @param coords 其他坐标（y:2D or y,z:3D）
	 */
	public Node(int globalIndex, double x, double ...coords) {
		this.globalIndex = globalIndex;
		this.coords[0] = x;
		if(coords!=null && coords.length > 0) {
			this.dim = 1+coords.length;
			for(int i=0;i<coords.length;i++)
				this.coords[1+i] = coords[i];
		} else {
			this.dim = 1;
		}
	}
	
	/**
	 * Set node globalIndex and coordinates
	 * 
	 * @param globalIndex
	 * @param coords
	 * @return
	 */
	public Node set(int globalIndex, double ...coords) {
		this.globalIndex = globalIndex;
		if(coords!=null && coords.length > 0) {
			this.dim = coords.length;
			for(int i=0;i<dim;i++)
				this.coords[i] = coords[i];
		}
		return this;
	}
	
	@Override
	public int getIndex() {
		return this.globalIndex;
	}
	
	@Override
	public int dim() {
		return dim;
	}
	
	@Override
	public double coord(int index) {
		return coords[index-1];
	}
	
	@Override
	public double[] coords() {
		double[] rlt;
		if(this.dim < 3) {
			rlt = new double[dim];
			for(int i=0;i<dim;i++)
				rlt[i] = this.coords[i];
		} else
			rlt = this.coords;
		return rlt;
	}
	
	public void setCoord(int index,double val) {
		coords[index-1] = val;
	}
	
	@Override
	public boolean coordEquals(Point p) {
		if(this.dim != p.dim())
			return false;
		for(int i=1;i<=this.dim;i++) {
			if(Math.abs(this.coord(i)-p.coord(i)) > Constant.meshEps)
				return false;
		}
		return true;
	}
	
	/**
	 * <P>Return true if the node is a inner node</P>
	 * <P>判断是否是内部结点</P>
	 * 
	 * <B>Note:</B>
	 * Function <code>Mesh.computeNodeBelongsToElements()</code> or 
	 * <code>Mesh.markBorderNode()</code> must be called before calling
	 * this function
	 * 
	 * @return
	 */
	public boolean isInnerNode() {
		if(this.getNodeType() == null) {//没有设置结点类型，按照几何结构判断
			
			if(belongToElements==null || belongToElements.size()==0)
				throw new FutureyeException("Call Mesh.computeNodeBelongsToElements() first!");
			double sum = 0.0;
			double coef = 0;
			if(dim==1) { //1D 按照所属单元个数判断
				if(belongToElements.size()==2)
					return true;
				else
					return false;
			} else if(dim==2) {//2D 按照包围结点单元是否形成圆周角2*PI
				//计算以node为顶点，周围单元结点与之形成的夹角角度，如果为360度（2*PI），就说明是内点
				coef = 2;
				for(int j=1;j<=belongToElements.size();j++) {
					sum += belongToElements.at(j).getAngleInElement2D(this);
				}
			} else if(dim==3) {//3D 按照单位球面三角形面积计算
				coef = 4;
				for(int j=1;j<=belongToElements.size();j++) {
					sum += belongToElements.at(j).getUnitSphereTriangleArea(this);
				}
			}
			if(Math.abs(sum-coef*Math.PI) <= Constant.meshEps)//2D=2*PI，3D=4*PI，是内部结点
				return true;
			else
				return false;
		} else {//Use the result of Mesh.markBorderNode() 按照结点类型标记判断
			
			for(int i=1;i<=this.nodeTypes.size();i++)
				if(this.nodeTypes.at(i) == null ||
						this.nodeTypes.at(i) != NodeType.Inner)
					return false;
			return true;
		}
	}
	
	/**
	 * Get node type for scalar valued function problems
	 * 
	 * @return
	 */
	public NodeType getNodeType() {
		if(nodeTypes == null) 
			return null;
		return nodeTypes.at(1);
	}
	
	/**
	 * Get node type at component <tt>nVVFComponent</tt> for vector valued function problems 
	 * 
	 * 对于向量值问题，每个函数分量在同一结点上的类型（例如速度各个分量、压强的边界类型）不一定相同，
	 * 该函数返回函数分量<tt>nVVFComponent</tt>对应的边界结点类型
	 * 
	 * @param nVVFComponent: Component of vector valued functions
	 * @return
	 */
	public NodeType getNodeType(int nVVFComponent) {
		if(nodeTypes == null) 
			return null;
		return nodeTypes.at(nVVFComponent);
	}
	
	/**
	 * Set node type for scalar valued function problems
	 * 
	 * @param nodeType: node type
	 */
	public void setNodeType(NodeType nodeType) {
		if(nodeTypes == null) 
			nodeTypes = new ObjVector<NodeType>();
		nodeTypes.setSize(1);
		nodeTypes.set(1,nodeType);
	}
	
	/**
	 * Set node type at component <tt>nVVFComponent<tt> for vector valued function problems
	 * 
	 * @param nVVFComponent: Component of vector valued function
	 * @param nodeType: Node type
	 */
	public void setNodeType(int nVVFComponent, NodeType nodeType) {
		if(nodeTypes == null) 
			nodeTypes = new ObjVector<NodeType>();
		if(nVVFComponent<=0) {
			throw new FutureyeException("nVVFComponent should be greater than 1!");
		}
		if(nodeTypes.size()<nVVFComponent)
			nodeTypes.setSize(nVVFComponent);
		nodeTypes.set(nVVFComponent,nodeType);
	}
	
	/**
	 * Add a element <code>e</code> which contains this node to 
	 * <code>this.belongToElements</code> list. 
	 * Duplicate element will NOT be added.
	 * 
	 * @param e
	 */
	public void addBelongToElements(Element e) {
		if(belongToElements == null)
			belongToElements = new ElementList();
		for(int i=1;i<=belongToElements.size();i++) {
			Element es = belongToElements.at(i);
			if(e.equals(es))//对象直接比较
				return;
			else if(e.globalIndex != 0 && e.globalIndex == es.globalIndex)//全局索引比较
				return;
		}
		belongToElements.add(e);
	}
	
	/**
	 * Add a <code>node</code> which is the neighbor of this node to 
	 * <code>this.neighbors</code> list.
	 * Duplicate node will NOT be added.
	 * 
	 * @param node
	 */
	public void addNeighbors(Node node) {
		if(neighbors == null)
			neighbors = new NodeList();
		for(int i=1;i<=neighbors.size();i++) {
			Node ni = neighbors.at(i);
			if(ni.equals(node))//对象直接比较
				return;
			else if(ni.globalIndex != 0 && ni.globalIndex == node.globalIndex)//全局索引比较
				return;
		}
		neighbors.add(node);		
	}
	
	public void clearBelongToElements() {
		if(belongToElements != null)
			belongToElements.clear();
	}
	
	public void clearNeighbors() {
		if(neighbors != null)
			neighbors.clear();
	}
	
	/**
	 * Get node level number from refined mesh
	 * 
	 * @return
	 */
	public int getRefineLevel() {
		return this.refineLevel;
	}
	
	/**
	 * Set node level number when doing mesh refinement
	 * 
	 * @param level
	 */
	public void setRefineLevel(int level) {
		this.refineLevel = level;
	}
	
	/**
	 * Return node information. Prefix "GN" means "Global Node"
	 * 
	 */
	public String toString()  {
		String s = "GN"+globalIndex+"( ";
		for(int i=0;i<dim;i++)
			s += String.valueOf(coords[i])+" ";
		return s+")";
	}
	
	@Override
    public boolean equals(Object obj) {
        if(super.equals(obj)) {
        	return true;
        } else {
        	if(obj instanceof Node) {
        		Node n2 = (Node)obj;
        		//if(globalIndex != 0 && globalIndex == n2.globalIndex && dim == n2.dim) {
        		if(globalIndex == n2.globalIndex && dim == n2.dim) {
        			//全局索引相同，坐标值相同
        			for(int i=0;i<dim;i++) {
        				if(Math.abs(this.coords[i] - n2.coords[i]) > Constant.meshEps) {
        					throw new FutureyeException(
        							this+" and "+n2+" should be equal, please check!");
        				}
        			}
        			return true;
        		}
        	}
        }
        return false;
    }
    
    @Override
    public int hashCode() {
    	if(globalIndex != 0) 
    		return globalIndex;
    	else 
    		return super.hashCode();
    }
}
 