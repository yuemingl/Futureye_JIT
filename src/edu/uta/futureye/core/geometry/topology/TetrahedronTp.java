package edu.uta.futureye.core.geometry.topology;

/**
 * 四面体拓扑结构
 * @author liuyueming
 *
 */
public class TetrahedronTp implements Topology3D {
	public static int[] vertices = {1,2,3,4};
	public static int[][] edges = {{1,2},{2,3},{3,1},{1,4},{2,4},{3,4}};
	public static int[][] faces = {{1,2,3},{1,2,4},{2,3,4},{3,1,4}};
	
	@Override
	public int[] getVertices() {
		return vertices;
	}
	
	@Override
	public int[][] getEdges() {
		return edges;
	}
	
	@Override
	public int[][] getFaces() {
		return faces;
	}
	
}
