package edu.uta.futureye.core;

import edu.uta.futureye.core.geometry.Point;
import edu.uta.futureye.util.Constant;
import edu.uta.futureye.util.list.ElementList;
import edu.uta.futureye.util.list.NodeList;

/**
 * Global finite element node class
 * 有限元（全局）结点
 * 
 * @author liuyueming
 *
 */
public class Node implements Point {
	//空间维度
	protected int dim = 0;
	
	//空间坐标
	protected double[] coords = new double[3];
	
	//全局索引（编号）
	public int globalIndex = 0;
	
	//结点类型：边界结点 or 内部结点
	private NodeType nodeType = null;
	
	//全局来看，结点所属单元
	public ElementList belongToElements = new ElementList();
	
	//相邻结点
	public NodeList neighbors = new NodeList();
	
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
		for(int i=1;i<=this.dim;i++) {
			if(Math.abs(this.coord(i)-p.coord(i)) > Constant.eps)
				return false;
		}
		return true;
	}
	
	public boolean isInnerNode() {
		return nodeType==NodeType.Inner;
	}
	
	public NodeType getNodeType() {
		return nodeType;
	}
	
	public void setNodeType(NodeType nodeType) {
		this.nodeType = nodeType;
	}
	
	public void addBelongToElements(Element e) {
		for(int i=1;i<=belongToElements.size();i++) {
			//TODO e.globalIndex 用索引代替对象直接比较？
			if(e.equals(belongToElements.at(i)))
				return;
		}
		belongToElements.add(e);
	}
	
	//////////////////////////////////////////////////////
	//结点属于的加密层次，即第level次加密产生的结点
	protected int level = 1;
	
	public int getLevel() {
		return this.level;
	}
	
	public void setLevel(int level) {
		this.level = level;
	}
	///////////////////////////////////////////////////////
	
	public String toString()  {
		String r = "GN"+globalIndex+"( ";
		for(int i=0;i<dim;i++)
			r += String.valueOf(coords[i])+" ";
		return r+")";
	}
}
 