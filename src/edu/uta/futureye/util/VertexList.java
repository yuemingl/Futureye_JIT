package edu.uta.futureye.util;

import java.util.LinkedList;
import java.util.List;

import edu.uta.futureye.core.Vertex;

public class VertexList {
	protected List<Vertex> vertices = new LinkedList<Vertex>();

	/**
	 * @param index start from 1,2,3...
	 * @return
	 */
	public Vertex at(int index) {
		if(index < 1)
			System.out.println("ERROR: VertexList index="+index);		
		return vertices.get(index-1);
	}

	public void add(Vertex edge) {
		this.vertices.add(edge);
	}
	
	public int size() {
		return vertices.size();
	}
	
	public void clear() {
		vertices.clear();
	}
	
	public String toString() {
		return vertices.toString();
	}	
}
