package edu.uta.futureye.util;

import java.util.ArrayList;
import java.util.List;

import edu.uta.futureye.core.Edge;

public class EdgeList {

	protected List<Edge> edges = new ArrayList<Edge>();

	/**
	 * @param index start from 1,2,3...
	 * @return
	 */
	public Edge at(int index) {
		if(index < 1)
			System.out.println("ERROR: EdgeList index="+index);		
		return edges.get(index-1);
	}

	public void add(Edge edge) {
		this.edges.add(edge);
	}
	
	public int size() {
		return edges.size();
	}
	
	public void clear() {
		edges.clear();
	}
	
	public String toString() {
		return edges.toString();
	}	
}

