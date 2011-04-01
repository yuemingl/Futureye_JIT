package edu.uta.futureye.core.geometry.topology;

/**
 * 三角形拓扑结构
 * @author liuyueming
 *
 */
public class TriangleTp implements Topology2D {
	public static int[] vertices = {1,2,3};
	public static int[][] edges = {{1,2},{2,3},{3,1}};
	//二次元
	//public static int[][] edges3 = {{1,2,4},{2,3,5},{3,1,6}};
	
	@Override
	public int[] getVertices() {
		return vertices;
	}
	
	@Override
	public int[][] getEdges() {
		return edges;
	}
}
