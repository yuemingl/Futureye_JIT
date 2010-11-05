package edu.uta.futureye.core;

import java.util.List;
import java.util.LinkedList;

import edu.uta.futureye.core.NodeType;
import edu.uta.futureye.core.intf.GeoEntity;
import edu.uta.futureye.core.intf.Point;
import edu.uta.futureye.util.DOFList;
import edu.uta.futureye.util.NodeList;
import edu.uta.futureye.util.Utils;

public class Edge implements GeoEntity{
	public Element owner = null;
	public List<NodeLocal> nodeLocalList = new LinkedList<NodeLocal>();
	
	public NodeType getBorderNodeType() {
		NodeType nt1 = nodeLocalList.get(0).globalNode.getNodeType();
		NodeType nt2 = nodeLocalList.get(nodeLocalList.size()-1).globalNode.getNodeType();
		
		if(nt1 == nt2) return nt1;
		else
			return null;
	}
	
	public boolean isBorderEdge() {
		for(NodeType nt : NodeType.values()) {
			if(nt == NodeType.Inner) continue;
			if(nt == this.getBorderNodeType()) return true;
		}
		return false;
	}
	
	/**
	 * Edge自己变为一个单元，用于边界积分（线积分）
	 * @return
	 */
	public Element changeToElement() {
		Element be = new Element();
		
		Node n1 = nodeLocalList.get(0).globalNode;
		Node n2 = nodeLocalList.get(1).globalNode;
		int ownerLocalIndex1 = this.owner.getLocalIndex(n1);
		int ownerLocalIndex2 = this.owner.getLocalIndex(n2);
		
		be.stemGeoEntity = this;
		
		be.addNode(n1,true);
		be.addNode(n2,true);
		
		
		DOFList eDOFList = owner.getDOFList(ownerLocalIndex1);
		for(int j=1;eDOFList!=null && j<=eDOFList.size();j++) {
			DOF dof = new DOF(
						1,
						n1.globalIndex,
						eDOFList.at(j).shapeFunction.restrictTo(1)
					);
			be.addDOF(1, dof);
		}
		eDOFList = owner.getDOFList(ownerLocalIndex2);
		for(int j=1;eDOFList!=null && j<=eDOFList.size();j++) {
			DOF dof = new DOF(
						2,
						n2.globalIndex,
						eDOFList.at(j).shapeFunction.restrictTo(2)
					);
			be.addDOF(2, dof);
		}	
		
		int dofIndex = 3;
		for(int i=2;i<nodeLocalList.size();i++) {
			NodeLocal nl = nodeLocalList.get(i);
			Node nn = nl.globalNode;
			be.addNode(nn,false);
			int ownerLocalIndexNN = this.owner.getLocalIndex(nn);
			eDOFList = owner.getDOFList(ownerLocalIndexNN);
			for(int j=1;eDOFList!=null && j<=eDOFList.size();j++) {
				int localIndex = dofIndex++;
				DOF dof = new DOF(
					localIndex,
					eDOFList.at(j).globalIndex,
					eDOFList.at(j).shapeFunction.restrictTo(localIndex)
				);
				be.addDOF(localIndex, dof);
			}
		}
		return be;
	}
	
	/**
	 * nodeLocalList 的前两个元素是Edge的端点
	 * @param coord
	 * @return
	 */
	public boolean isCoordOnEdge(double coord[]) {
		Node n1 = nodeLocalList.get(0).globalNode;
		Node n2 = nodeLocalList.get(1).globalNode;
		Node n = new Node(coord);
		return Utils.isPointOnLineSegment(n1,n2,n);
	}
	
	public NodeList getEndNodes() {
		NodeList rlt = new NodeList();
		rlt.add(nodeLocalList.get(0).globalNode);
		rlt.add(nodeLocalList.get(1).globalNode);
		return rlt;
	}
	
	public String toString() {
		return nodeLocalList.toString();
	}
}
