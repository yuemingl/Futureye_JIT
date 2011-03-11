package edu.uta.futureye.core;

import edu.uta.futureye.core.geometry.Point;

/**
 * * Local finite element node class. 
 * 
 * * 有限元（局部）结点，包含对全局结点的引用，其坐标值为全局结点的坐标。
 * * 局部结点的内存占用相对全局结点要少很多
 * * 局部结点的局部索引（编号）是相对于该节点所属的对象来说的，例如：
 *     三角形单元，其上的局部结点为t1、t2和t3，有局部边e12、e23和e31，
 *     考察e23对应的全局边ge23，其上的局部结点 编号应该是ge1、ge2，
 *     但他们相关联的全局结点有关系：t2=ge1，t3=ge2
 * * 当一个局部结点属于某个有限单元时，局部结点应该包含对单元的引用，
 *   不属于某个有限单元时，该应用可以为null，例如一个全局边上的局部结点
 *   就不属于任何有限单元
 * 
 * @author liuyueming
 *
 */
public class NodeLocal implements Point {
	public Element owner = null; //所属单元
	public Node globalNode; //对应相同坐标的全局结点
	public int localIndex; //在单元内的局部编号
	
	public NodeLocal(int localIndex, Node globalNode) {
		this.localIndex = localIndex;
		this.globalNode = globalNode;
	}
	
	public NodeLocal(int localIndex, Node globalNode, Element owner) {
		this.localIndex = localIndex;
		this.globalNode = globalNode;
		this.owner = owner;
	}
	
	public Node globalNode() {
		return this.globalNode;
	}
	
	@Override
	public double coord(int index) {
		return this.globalNode.coord(index);
	}

	@Override
	public boolean coordEquals(Point p) {
		return this.globalNode.coordEquals(p);
	}

	@Override
	public double[] coords() {
		return this.globalNode.coords();
	}

	@Override
	public int dim() {
		return this.globalNode.dim();
	}

	@Override
	public int getIndex() {
		return this.localIndex;
	}

	@Override
	public void setCoord(int index, double value) {
		this.globalNode.setCoord(index, value);
	}
	
	@Override
	public String toString() {
		return "NodeLocal"+localIndex+"<=>"+globalNode.toString();
	}
}
