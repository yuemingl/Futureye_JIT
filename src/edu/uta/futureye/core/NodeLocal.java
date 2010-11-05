package edu.uta.futureye.core;

public class NodeLocal {
	
	public NodeLocal(Node globalNode, int localIndex) {
		this.globalNode = globalNode;
		this.localIndex = localIndex;
	}
	
	public Node globalNode;
	public int localIndex;
	
	public String toString() {
		return "LN"+localIndex+":"+globalNode.toString();
	}
}
