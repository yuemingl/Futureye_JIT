package edu.uta.futureye.core;

import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.core.geometry.GeoEntity2D;
import edu.uta.futureye.util.FutureyeException;
import edu.uta.futureye.util.container.DOFList;
import edu.uta.futureye.util.container.ObjList;
import edu.uta.futureye.util.container.VertexList;

/**
 * Local face of an element
 * 局部面（三维单元的局部面）
 * 
 * @author liuyueming
 *
 */
public class FaceLocal extends GeoEntity2D<EdgeLocal,NodeLocal> {
	public int localIndex;
	//global face shared with all elements that containing the face
	protected Face globalFace = null; 
	public Element owner = null;

	protected Vector localUnitNormVector;

	public Vector getLocalUnitNormVector() {
		return localUnitNormVector;
	}

	public void setLocalUnitNormVector(Vector localUnitNormVector) {
		this.localUnitNormVector = localUnitNormVector;
	}

	public FaceLocal(int localIndex, Element owner) {
		this.localIndex = localIndex;
		this.owner = owner;
	}
	
	public FaceLocal(int localIndex, Face globalFace) {
		this.localIndex = localIndex;
		this.globalFace = globalFace;
	}
	
	public FaceLocal(int localIndex, Face globalFace, Element owner) {
		this.localIndex = localIndex;
		this.globalFace = globalFace;
		this.owner = owner;
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
	
	public Face getGlobalFace() {
		return globalFace;
	}

	public void setGlobalFace(Face globalFace) {
		this.globalFace = globalFace;
	}
	
	public Face buildFace() {
		Face face = new Face();
		ObjList<Vertex> vertices = this.getVertices();
		for(int n=1;n<=vertices.size();n++) {
			face.addVertex(new Vertex(n,new NodeLocal(n,vertices.at(n).globalNode())));
		}
		ObjList<EdgeLocal> localEdges = this.getEdges();
		for(int n=1;n<=localEdges.size();n++) {
			face.addEdge(new EdgeLocal(n, localEdges.at(n).globalEdge));
		}
		ObjList<NodeLocal> localFaceNodes = this.getFaceNodes();
		if(localFaceNodes != null && localFaceNodes.size()>0) {
			for(int n=1;n<=localFaceNodes.size();n++) {
				face.addFaceNode(new NodeLocal(n,localFaceNodes.at(n).globalNode));
			}
		}
		face.setTopology(this.getTopology());
		return face;
	}
	
	/**
	 * Face自己变成二维单元，用于边界积分（面积分）
	 * @return
	 */
	public Element changeToElement() {
		
		Element be = new Element(this.buildFace());
		
		//Node DOFs
		VertexList beVertices = be.vertices();
		int localDOFIndex = 1;
		for(int i=1;i<=beVertices.size();i++) {
			DOFList eDOFList = owner.getNodeDOFList(beVertices.at(i).localNode().localIndex);
			for(int j=1;eDOFList!=null && j<=eDOFList.size();j++) {
				DOF dof = new DOF(
							localDOFIndex,
							eDOFList.at(j).globalIndex,
							eDOFList.at(j).getSSF().restrictTo(localDOFIndex)
						);
				be.addNodeDOF(localDOFIndex, dof);
				localDOFIndex++;
			}
			
		}
		//TODO Edge DOFs?
		//TODO Face DOFs?
		
		return be;		
		
	}
	
	public String toString() {
		if(this.globalFace != null)
			return "LocalFace"+localIndex+"<=>"+globalFace.toString();
		else
			return "LocalFace"+localIndex+this.vertices.toString();
	}
}