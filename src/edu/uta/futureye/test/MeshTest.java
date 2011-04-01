package edu.uta.futureye.test;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.io.MeshReader;
import edu.uta.futureye.util.container.EdgeList;
import edu.uta.futureye.util.container.ElementList;
import edu.uta.futureye.util.container.NodeList;

public class MeshTest {
	public static void main(String[] args) {
		test1();
		test2();
	}
	
	public static void test1() {
	    MeshReader reader = new MeshReader("patch_triangle.grd");
	    Mesh mesh = reader.read2DMesh();		
	    NodeList nodes = mesh.getNodeList();
	    ElementList elements = mesh.getElementList();
	    
	    for(int i=1;i<=nodes.size();i++) {
	    	System.out.println(nodes.at(i));
	    }
	    for(int i=1;i<=elements.size();i++) {
	    	System.out.println(elements.at(i));
	    }
	    
	    Node node = nodes.at(9);
	    System.out.println("GN"+node.globalIndex+" coord:"+
	    		node.coord(1)+","+node.coord(2));
	    double []coord = node.coords();
	    System.out.println("GN"+node.getIndex()+" coord:"+
	    		coord[0]+","+coord[1]);
	    
	    Element ele = elements.at(1);
	    System.out.println("GE"+ele.globalIndex+" nodes:"+
	    		ele.nodes);
	    
	    List<Node> list = nodes.toList();
	    Collections.sort(list, new Comparator<Node>(){
	        @Override
	        public int compare(Node o1, Node o2) {
	        	//desc
	            if(o1.globalIndex > o2.globalIndex)
	                return -1;
	            else
	                return 1;
	        }
	    });
	    for(int i=1;i<=nodes.size();i++) {
	    	System.out.println(nodes.at(i));
	    }
	    
	}
	
	public static void test2() {
	    MeshReader reader = new MeshReader("patch_triangle.grd");
	    Mesh mesh = reader.read2DMesh();		
	    NodeList nodes = mesh.getNodeList();
	    ElementList elements = mesh.getElementList();

	    mesh.computeNodeBelongsToElements();
	    //Print all elements that a node is contained in them
	    for(int i=1;i<=nodes.size();i++) {
	    	System.out.println(nodes.at(i)+": "+
	    			nodes.at(i).belongToElements);
	    }
	    
	    mesh.computeNeighborNodes();
	    //Print all neighbors of a node
	    for(int i=1;i<=nodes.size();i++) {
	    	System.out.println(nodes.at(i)+": "+
	    			nodes.at(i).neighbors);
	    }
	    
	    mesh.computeGlobalEdge();
	    EdgeList edges = mesh.getEdgeList();
	    mesh.computeNeighborElements();
	    //Print all global edges
	    for(int i=1;i<=edges.size();i++) {
	    	System.out.println(edges.at(i));
	    }
	    //Print all neighbors of an element
	    for(int i=1;i<=elements.size();i++) {
	    	System.out.println(elements.at(i)+": "+
	    			elements.at(i).neighbors);
	    }
	}
}
