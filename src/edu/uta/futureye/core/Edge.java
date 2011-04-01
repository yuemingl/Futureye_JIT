package edu.uta.futureye.core;

import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.core.geometry.GeoEntity1D;
import edu.uta.futureye.util.Utils;
import edu.uta.futureye.util.container.NodeList;

/**
 * Global edge of an element
 * 全局边
 * 
 * @author liuyueming
 *
 */
public class Edge extends GeoEntity1D<NodeLocal> {
	protected int globalIndex;
	private Vector globalUnitNormVector; //全局单位法方向
	
	public Edge() {
	}
	
	public Edge(NodeLocal begin, NodeLocal end) {
		addVertex(new Vertex(1,begin));
		addVertex(new Vertex(2,end));

		globalUnitNormVector = Utils.getNormVector(
				this.vertices.at(1),
				this.vertices.at(2)
				);
	}
	
	public int getGlobalIndex() {
		return globalIndex;
	}

	public void setGlobalIndex(int globalIndex) {
		this.globalIndex = globalIndex;
	}
	
	public Vector getNormVector() {
		if(this.globalUnitNormVector == null)
			globalUnitNormVector = Utils.getNormVector(
					this.vertices.at(1),
					this.vertices.at(2)
					);
			
		return this.globalUnitNormVector;
	}
	
    public NodeType getBorderType() {                 
    	return getBorderType(1);
    }
    
	/**
	 * 对于向量值问题，每个分量在同一边界上的类型不一定相同，
	 * 该函数返回分量<tt>vvfIndex</tt>对应的边界类型
	 * Vector valued function (vvf)
	 * @param vvfIndex
	 * @return
	 */
    public NodeType getBorderType(int vvfIndex) {                 
    	NodeType nt1 = this.vertices.at(1).globalNode().getNodeType(vvfIndex);
    	NodeType nt2 = this.vertices.at(2).globalNode().getNodeType(vvfIndex);                       
    	if(nt1 == nt2) return nt1;                  
    	else {
    		//TODO Exception?
    		return null;
    	}
    }
    
	public boolean isBorderEdge() {
		if(this.vertices.at(1).globalNode().getNodeType()==NodeType.Inner ||
				this.vertices.at(2).globalNode().getNodeType()==NodeType.Inner)
			return false;
		else
			return true;
		
//		for(NodeType nt : NodeType.values()) {
//			if(nt == NodeType.Inner) continue;
//			if(nt == this.edgeNodes.at(1).getNodeType() &&
//					nt == this.edgeNodes.at(this.edgeNodes.size()).getNodeType()) 
//				return true;
//		}
//		return false;
	}
	
	/**
	 * 判断某个坐标点是否在edge上
	 * @param coord
	 * @return
	 */
	public boolean isCoordOnEdge(double coords[]) {
		return Utils.isPointOnLineSegment(
				this.vertices.at(1),
				this.vertices.at(2),
				new Vertex().set(0, coords));
	}
	
    public Vertex beginVertex() {
		return this.vertices.at(1);
    }
    
    public Vertex endVertex() {
		return this.vertices.at(2);
    }
    
	public Node beginNode() {
		return this.vertices.at(1).globalNode();
	}

	public Node endNode() {
		return this.vertices.at(2).globalNode();
	}
	
	public NodeList getEndNodes() {
		NodeList rlt = new NodeList();
		rlt.add(this.vertices.at(1).globalNode());
		rlt.add(this.vertices.at(2).globalNode());
		return rlt;
	}
	
	public double getEdgeLength() {
		return Utils.computeLength(
					this.vertices.at(1),
					this.vertices.at(2)
				);
	}
	
	public String toString() {
		return "Edge"+this.globalIndex+":"+this.vertices.toString();
	}	
}
