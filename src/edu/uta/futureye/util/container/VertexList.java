package edu.uta.futureye.util.container;

import edu.uta.futureye.core.Vertex;

/**
 * Vertex List Class
 * 顶点类列表
 * 
 * @author liuyueming
 *
 */
public class VertexList extends ObjList<Vertex> {
	@Override
	public String toString() {
		return "VertexList"+objs.toString();
	}
	
}
