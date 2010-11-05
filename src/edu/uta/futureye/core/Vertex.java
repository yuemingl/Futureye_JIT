package edu.uta.futureye.core;

import edu.uta.futureye.core.intf.Point;
import edu.uta.futureye.util.Constant;

/**
 * 单元顶点
 * @author liuyueming
 *
 */
public class Vertex implements Point{
	public Element owner = null;
	int dim = 0;
	protected double[] coords = new double[3];
	public int localIndex;
	
	public Vertex(int dim) {
		this.dim = dim;	
	}
	
	public void set(int localIndex, double ...coords) {
		this.localIndex = localIndex;
		this.coords = coords;
	}
	
	/**
	 * 将一个Node坐标赋值给Vertex
	 * @param localIndex
	 * @param node
	 */
	public void set(int localIndex, Node node) {
		this.localIndex = localIndex;
		for(int i=0;i<node.dim();i++)
			coords[i] = node.coord(i+1);
	}

	@Override
	public double coord(int index) {
		return coords[index-1];
	}

	@Override
	public int dim() {
		return dim;
	}
	
	public String toString() {
		String s = "";
		if(this.owner != null)
			s += "GN"+owner.nodes.at(localIndex).globalIndex;
		s += "( ";
		for(int i=0;i<dim;i++)
			s += String.valueOf(coords[i])+" ";
		return s+")";
	}

	public boolean coordEquals(Point p) {
		for(int i=1;i<=this.dim;i++) {
			if(Math.abs(this.coord(i)-p.coord(i)) > Constant.eps)
				return false;
		}
		return true;
	}
}
