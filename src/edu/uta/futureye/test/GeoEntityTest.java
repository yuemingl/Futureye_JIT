package edu.uta.futureye.test;

import edu.uta.futureye.core.EdgeLocal;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.core.NodeLocal;
import edu.uta.futureye.core.Vertex;
import edu.uta.futureye.core.geometry.GeoEntity2D;
import edu.uta.futureye.util.container.NodeList;

public class GeoEntityTest {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		NodeList nodes = new NodeList();
		nodes.add(new Node(1, 0.0,1.0));
		nodes.add(new Node(2, 0.0,0.0));
		nodes.add(new Node(3, 1.0,0.0));
		nodes.add(new Node(4, 1.0,1.0));
		nodes.add(new Node(5, 0.0,0.5));
		nodes.add(new Node(6, 0.5,0.0));
		nodes.add(new Node(7, 1.0,0.5));
		nodes.add(new Node(8, 0.5,1.0));
		nodes.add(new Node(9, 0.5,0.5));

		Element e = new Element();
		GeoEntity2D<EdgeLocal,NodeLocal> face = 
			new GeoEntity2D<EdgeLocal,NodeLocal>();
		//vertices of face: 1,2,3,4
		for(int i=1;i<=4;i++)
			face.addVertex(new Vertex(i,new NodeLocal(i,nodes.at(i))));
		
		int[] idxLoop = {1,2,3,4,1};
		int[] edgeIdx = {5,6,7,8};
		for(int i=0;i<4;i++) {
			EdgeLocal el = new EdgeLocal(i,e);
			int idxBeg = idxLoop[i];
			int idxEnd = idxLoop[i+1];
			//vertices of edge: idxBeg,idxEnd
			el.addVertex(new Vertex(idxBeg,new NodeLocal(idxBeg,nodes.at(idxBeg))));
			el.addVertex(new Vertex(idxEnd,new NodeLocal(idxEnd,nodes.at(idxEnd))));
			//node on edge
			el.addEdgeNode(new NodeLocal(edgeIdx[i],nodes.at(edgeIdx[i])));
			//edge on face
			face.addEdge(el);
		}
		//node on face
		face.addFaceNode(new NodeLocal(9,nodes.at(9)));
		//create the element
		e.setGeoEntity(face);
		//output info
		System.out.println(e);
		for(int i=1;i<=e.nodes.size();i++)
			System.out.println(e.nodes.at(i));
	}

}
