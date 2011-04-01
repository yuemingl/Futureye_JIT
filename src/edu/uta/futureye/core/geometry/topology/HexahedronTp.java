package edu.uta.futureye.core.geometry.topology;

/**
 * 六面体拓扑结构
 * @author liuyueming
 *
 */
public class HexahedronTp implements Topology3D {
	public static int[] vertices = {1,2,3,4,5,6,7,8};
	public static int[][] edges = {
		{1,2},{2,3},{3,4},{4,1},
		{5,6},{6,7},{7,8},{8,5},
		{1,5},{2,6},{3,7},{4,8}};
	public static int[][] faces = {
		{1,2,3,4},{8,7,6,5},
		{2,1,5,6},{3,2,6,7},{4,3,7,8},{1,4,8,5}};
	
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
