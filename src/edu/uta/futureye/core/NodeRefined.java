package edu.uta.futureye.core;

import edu.uta.futureye.util.NodeList;

public class NodeRefined extends Node {
	public NodeList constrainNodes = new NodeList();
	
	public NodeRefined(int dim) {
		super(dim);
		this.level = 2;
	}
	
	public void addConstrainNode(Node node) {
		for(int i=1;i<=this.constrainNodes.size();i++)
			if(node.equals(this.constrainNodes.at(i)))
				return;
		constrainNodes.add(node);
	}
	
	public void clearConstrainNode() {
		this.constrainNodes.clear();
	}	
	
	/**
	 * 判断是否为Hanging node
	 * @return
	 */
	public boolean isHangingNode() {
		//没有constrain node的加密结点不是hanging node
		return this.constrainNodes.size()>0;
	}
	

}
