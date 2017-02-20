/**
 * Copyright (c) 2010, nkliuyueming@gmail.com. All rights reserved.
 * 
 * 
 */
package edu.uta.futureye.core;

import edu.uta.futureye.algebra.SpaceVector;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.core.geometry.GeoEntity2D;
import edu.uta.futureye.core.geometry.Point;
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
			ObjList<Vertex> vs = this.getVertices();
			Point p1 = vs.at(1);
			Point p2 = vs.at(2);
			Point p3 = vs.at(3);
			SpaceVector s1 = new SpaceVector(3);
			SpaceVector s2 = new SpaceVector(3);
			s1.set(1, p2.coord(1)-p1.coord(1));
			s1.set(2, p2.coord(2)-p1.coord(2));
			s1.set(3, p2.coord(3)-p1.coord(3));
			
			s2.set(1, p3.coord(1)-p2.coord(1));
			s2.set(2, p3.coord(2)-p2.coord(2));
			s2.set(3, p3.coord(3)-p2.coord(3));
			
			Vector rlt = s1.crossProduct(s2);
			rlt.scale(1.0/rlt.norm2());
			return rlt;
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
	 * For vector valued problems, return boundary type of component <tt>nVVFComponent</tt>
	 * <p>
	 * 对于向量值问题，每个分量在同一边界上的类型不一定相同，
	 * 该函数返回分量<tt>nVVFComponent</tt>对应的边界类型
	 * 
	 * @param nVVFComponent
	 * @return
	 */
	public NodeType getBorderType(int nVVFComponent) {
		NodeType nt1 = this.vertices.at(1).globalNode()
				.getNodeType(nVVFComponent);
		for (int i = 2; i < this.vertices.size(); i++) {
			NodeType nt2 = this.vertices.at(2).globalNode()
					.getNodeType(nVVFComponent);
			if (nt1 != nt2)
				return null;
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
