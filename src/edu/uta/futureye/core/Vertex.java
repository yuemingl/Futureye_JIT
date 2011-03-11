package edu.uta.futureye.core;

import edu.uta.futureye.core.geometry.Point;
import edu.uta.futureye.util.Constant;

/**
 * 单元顶点
 * @author liuyueming
 *
 */
public class Vertex implements Point {
	protected int dim = 0;
	protected double[] coords = new double[3];

	protected int localIndex;
	public Element owner = null;
	//顶点对应的结点
	protected NodeLocal refNodeLocal = null;
	
	public Node globalNode() {
		return refNodeLocal.globalNode;
	}
	
	public NodeLocal localNode() {
		return refNodeLocal;
	}

	public Vertex() {
		
	}
	
	/**
	 * 构造一个顶点
	 * @param localIndex 局部索引（编号）
	 * @param x 第一个坐标
	 * @param coords 其他坐标（y:2D or y,z:3D）
	 */
	public Vertex(int localIndex, double x,double ...coords) {
		this.localIndex = localIndex;
		this.coords[0] = x;
		if(coords!=null && coords.length > 0) {
			this.dim = 1+coords.length;
			for(int i=0;i<dim;i++)
				this.coords[1+i] = coords[i];
		} else {
			this.dim = 1;
		}
	}
	
	/**
	 * 从一个NodeLocal构造Vertex，并与NodeLocal关联
	 * @param localIndex
	 * @param localNode
	 */
	public Vertex(int localIndex, NodeLocal localNode) {
		set(localIndex,localNode);
	}
	
	/**
	 * 置顶点编号和坐标值
	 * @param localIndex
	 * @param coords
	 * @return self
	 */
	public Vertex set(int localIndex, double ...coords) {
		this.localIndex = localIndex;
		if(coords!=null && coords.length > 0) {
			this.dim = coords.length;
			for(int i=0;i<dim;i++)
				this.coords[i] = coords[i];
		}
		return this;
	}
	
	/**
	 * 将一个NodeLocal坐标赋值给Vertex，并与NodeLocal关联
	 * @param localIndex
	 * @param node
	 * @return self
	 */
	public Vertex set(int localIndex, NodeLocal localNode) {
		this.localIndex = localIndex;
		this.dim = localNode.dim();
		for(int i=0;i<dim;i++)
			coords[i] = localNode.coord(i+1);
		this.refNodeLocal = localNode;
		return this;
	}
	
	/**
	 * get local index
	 */
	@Override
	public int getIndex() {
		return this.localIndex;
	}
	
	public void setLocalIndex(int index) {
		this.localIndex = index;
	}
	
	public void setNodeLocalIndex(int index) {
		this.refNodeLocal.localIndex = index;
	}
	
	public void setAllLocalIndex(int index) {
		this.localIndex = index;
		this.refNodeLocal.localIndex = index;
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
	
	@Override
	public void setCoord(int index, double value) {
		this.coords[index-1] = value;
	}
	
	@Override
	public boolean coordEquals(Point p) {
		for(int i=1;i<=this.dim;i++) {
			if(Math.abs(this.coord(i)-p.coord(i)) > Constant.eps)
				return false;
		}
		return true;
	}
	
	public String toString() {
		String s = "V"+this.localIndex;
		if(this.refNodeLocal != null)
			s += "<=>GN"+refNodeLocal.globalNode.globalIndex;
		s += "( ";
		for(int i=0;i<dim;i++)
			s += String.valueOf(coords[i])+" ";
		return s+")";
	}
}
