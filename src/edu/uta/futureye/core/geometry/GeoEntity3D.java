package edu.uta.futureye.core.geometry;

import edu.uta.futureye.core.geometry.topology.Topology3D;
import edu.uta.futureye.util.list.ObjList;

/**
 * 三维几何实体，保存有限元Element的几何信息
 * 
 *   父类的vertices：三维单元对应顶点集合（四面体单元：四个顶点，六面体单元：八个顶点，...）
 *   faces：三维单元的面集合
 *   volumeNodes：三维单元内部的结点集合
 *   
 * @author liuyueming
 *
 */
public class GeoEntity3D<
						TFace extends GeoEntity2D<TEdge,TNode>,
						TEdge extends GeoEntity1D<TNode>,
						TNode extends Point
						> extends GeoEntity0D {
	protected Topology3D topology = null;
	
	protected ObjList<TFace> faces = new ObjList<TFace>();
	protected ObjList<TNode> volumeNodes = null;
	
	public void addFace(TFace edge) {
		this.faces.add(edge);
	}
	public void addAllFaces(ObjList<TFace> edges) {
		this.faces.clear();
		this.faces.addAll(edges);
	}
	public ObjList<TFace> getFaces() {
		return this.faces;
	}
	public void clearFaces() {
		this.faces.clear();
	}
	
	public void addVolumeNode(TNode node) {
		if(this.volumeNodes == null)
			this.volumeNodes = new ObjList<TNode>();
		this.volumeNodes.add(node);
	}
	public void addAllVolumeNodes(ObjList<TNode> faceNodes) {
		if(this.volumeNodes == null)
			this.volumeNodes = new ObjList<TNode>();
		this.volumeNodes.clear();
		this.volumeNodes.addAll(faceNodes);
	}
	public ObjList<TNode> getVolumeNodes() {
		return this.volumeNodes;
	}
	public void clearVolumeNodes() {
		if(this.volumeNodes != null)
			this.volumeNodes.clear();
	}
	public void clearAll() {
		this.faces.clear();
		if(this.volumeNodes != null)
			this.volumeNodes.clear();
	}
	public Topology3D getTopology() {
		return topology;
	}
	public void setTopology(Topology3D topology) {
		this.topology = topology;
	}		
}
