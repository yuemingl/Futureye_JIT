package edu.uta.futureye.core.geometry;

import edu.uta.futureye.util.list.ObjList;

/**
 * 一维几何实体，保存有限元Element的几何信息
 *   父类的vertices：一维单元对应的线段顶点（两个端点）
 *   nodes：边上的结点列表，不包括端点（e.g.线性元：0个，二次元：1个，...）
 *   
 * @author liuyueming
 *
 */
public class GeoEntity1D<TNode extends Point> extends GeoEntity0D {
	//边上的结点，不包括端点
	protected ObjList<TNode> edgeNodes = null;
	
	public void addEdgeNode(TNode node) {
		if(this.edgeNodes == null)
			this.edgeNodes = new ObjList<TNode>();
		this.edgeNodes.add(node);
	}
	
	public void addAllEdgeNodes(ObjList<TNode> nodes) {
		if(this.edgeNodes == null)
			this.edgeNodes = new ObjList<TNode>();
		this.edgeNodes.clear();
		this.edgeNodes.addAll(nodes);
	}
	
	public ObjList<TNode> getEdgeNodes() {
		return this.edgeNodes;
	}
	
	public void clearEdgeNodes() {
		if(this.edgeNodes != null)
			this.edgeNodes.clear();
	}
}