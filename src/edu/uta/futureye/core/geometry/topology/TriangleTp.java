package edu.uta.futureye.core.geometry.topology;

/**
 * 三角形拓扑结构
 * @author liuyueming
 *
 */
public class TriangleTp implements Topology2D {
	public static int[] vertices = {1,2,3};
	public static int[][] edges = {{1,2},{2,3},{3,1}};
	
	@Override
	public int[] getVertices() {
		return vertices;
	}
	
	@Override
	public int[][] getEdges() {
		return edges;
	}
}
