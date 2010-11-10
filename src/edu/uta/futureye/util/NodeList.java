package edu.uta.futureye.util;

import java.util.ArrayList;
import java.util.List;

import edu.uta.futureye.core.Node;

/**
 * Node List Class
 * 节点列表类
 * @author liuyueming
 *
 */
public class NodeList {
	protected List<Node> nodes = new ArrayList<Node>();

	/**
	 * @param index start from 1,2,3...
	 * @return
	 */
	public Node at(int index) {
		if(index < 1)
			System.out.println("ERROR: NodeList index="+index);		
		return nodes.get(index-1);
	}

	public void add(Node node) {
		this.nodes.add(node);
	}
	
	public void addAll(NodeList list) {
		this.nodes.addAll(list.nodes);
	}
	
	public int size() {
		return nodes.size();
	}
	
	public void clear() {
		nodes.clear();
	}
	
	public String toString() {
		return nodes.toString();
	}
	
}
