package edu.uta.futureye.core;

import edu.uta.futureye.core.intf.Point;
import edu.uta.futureye.util.Constant;
import edu.uta.futureye.util.ElementList;
import edu.uta.futureye.util.NodeList;

public class Node implements Point {
	public ElementList belongToElements = new ElementList();
	public NodeList neighbors = new NodeList();
	public int globalIndex = 0;
	
	private int dim = 0;
	private double[] coords = new double[3];
	
	private NodeType nodeType = null;
	
	public Node(int dim) {
		this.dim = dim;
	}

	/**
	 * 需要考虑该构造函数是否有必要
	 * @param coords
	 */
	public Node(double[] coords) {
		this.dim = coords.length;
		for(int i=0;i<this.dim;i++) {
			this.coords[i] = coords[i];
		}
	}
	
	public void set(int globalIndex, double ...coords) {
		this.globalIndex = globalIndex;
		this.coords = coords;
	}
	
	@Override
	public double coord(int index) {
		return coords[index-1];
	}
	
	public double[] coords() {
		return this.coords;
	}
	
	public void setCoord(int index,double val) {
		coords[index-1] = val;
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

	@Override
	public int dim() {
		return dim;
	}
	
	public void addBelongToElements(Element e) {
		for(int i=1;i<=belongToElements.size();i++) {
			//TODO ??? e.globalIndex ???
			if(e.equals(belongToElements.at(i)))
				return;
		}
		belongToElements.add(e);
	}
	
	public String toString()  {
		String r = "GN"+globalIndex+"( ";
		for(int i=0;i<dim;i++)
			r += String.valueOf(coords[i])+" ";
		return r+")";
	}
	
	public boolean coordEquals(Point p) {
		for(int i=1;i<=this.dim;i++) {
			if(Math.abs(this.coord(i)-p.coord(i)) > Constant.eps)
				return false;
		}
		return true;
	}
	
	//////////////////////////////////////////////////////
	//加密层次
	protected int level = 1;
	
	public int getLevel() {
		return this.level;
	}
	
	public void setLevel(int level) {
		this.level = level;
	}
	///////////////////////////////////////////////////////
}
 