package edu.uta.futureye.core;

import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.core.geometry.GeoEntity2D;
import edu.uta.futureye.util.FutureyeException;
import edu.uta.futureye.util.container.ObjList;

/**
 * Global face of an element
 * 全局面（三维单元的全局面）
 * 
 * @author liuyueming
 *
 */
public class Face extends GeoEntity2D<EdgeLocal,NodeLocal> {
	protected int globalIndex;
	protected Vector globalUnitNormVector = null; //全局单位法方向
	
	public int getGlobalIndex() {
		return globalIndex;
	}

	public void setGlobalIndex(int globalIndex) {
		this.globalIndex = globalIndex;
	}
	
	public Vector getNormVector() {
		if(this.globalUnitNormVector == null) {
			//TODO
		}
		return this.globalUnitNormVector;
	}
	
	/**
	 * 返回面的边界类型，确保所有顶点的类型要都相同
	 * 
	 * @return
	 */
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
    	for(int i=2;i<this.vertices.size();i++) {
    		NodeType nt2 = this.vertices.at(2).globalNode().getNodeType(vvfIndex);
    	   	if(nt1 != nt2) return null;                  
    	}
    	return nt1;                  
    }
    
	public boolean isBorderFace() {
		//顶点对应的NodeLocal是非Inner即可，也就是说只要有一个是Inner说明该面不是边界面
		ObjList<Vertex> vs = this.getVertices();
		if(vs.size() >= 3) {
			for(int i=1;i<=vs.size();i++) {
				NodeLocal nl = vs.at(i).localNode();
				if(nl.globalNode.getNodeType()==NodeType.Inner)
					return false;
			}
		} else {
			FutureyeException ex = new FutureyeException("Number of vertices on a face is "+vs.size());
			ex.printStackTrace();
			System.exit(0);
		}
		return true;
	}
	
	public String toString() {
		return "GlobalFace"+this.globalIndex+this.vertices.toString();
	}
}
