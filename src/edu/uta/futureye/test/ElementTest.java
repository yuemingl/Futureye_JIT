package edu.uta.futureye.test;

import edu.uta.futureye.core.Edge;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.util.EdgeList;

public class ElementTest {
	
	public static void testTriangleNode6() {
		/*
		 *  3
		 *  | \
		 *  |  \
		 *  6   5
		 *  |    \
		 *  |     \
		 *  1--4---2 
		 */
		Node[] nodes = new Node[6];
		for(int i=0;i<nodes.length;i++)
			nodes[i] = new Node(2);
		
		nodes[0].set(1, 0.0,0.0);
		nodes[1].set(2, 1.0,0.0);
		nodes[2].set(3, 0.0,1.0);
		nodes[3].set(4, 0.5,0.0);
		nodes[4].set(5, 0.5,0.5);
		nodes[5].set(6, 0.0,0.5);
		
		Element e = new Element();
		for(int i=0;i<nodes.length;i++) {
			e.addNode(nodes[i], i<3?true:false);
		}
		
		e.updateJacobinLinear2D();
		System.out.println(e.getJacobin());
		
		EdgeList eList = e.getEdgeList();
		for(int i=1;i<=eList.size();i++) {
			Edge edge = eList.at(i);
			System.out.println(edge);
		}
	}
	
	public static void testTriangleNode9() {
		/*
		 *  3
		 *  | \
		 *  |  \
		 *  6   8
		 *  |    \
		 *  9     5
		 *  |      \
		 *  |       \
		 *  1--4--7--2 
		 */
		Node[] nodes = new Node[9];
		for(int i=0;i<nodes.length;i++)
			nodes[i] = new Node(2);
		
		double onet = 1.0/3.0;
		double twot = 2.0/3.0;
		
		nodes[0].set(1, 0.0,0.0);
		nodes[1].set(2, 1.0,0.0);
		nodes[2].set(3, 0.0,1.0);
		nodes[3].set(4, onet,0.0);
		nodes[4].set(5, twot,onet);
		nodes[5].set(6, twot,0.0);
		nodes[6].set(7, twot,0.0);
		nodes[7].set(8, onet,twot);
		nodes[8].set(9, onet,0.0);
		
		Element e = new Element();
		for(int i=0;i<nodes.length;i++) {
			e.addNode(nodes[i], i<3?true:false);
		}
		
		EdgeList eList = e.getEdgeList();
		for(int i=1;i<=eList.size();i++) {
			Edge edge = eList.at(i);
			System.out.println(edge);
		}		
	}	
	
	public static void testRectangleNode8() {
		/*
		 * 4--7--3
		 * |     |
		 * 8     6
		 * |     |
		 * 1--5--2
		 * 
		 */	
		Node[] nodes = new Node[8];
		for(int i=0;i<nodes.length;i++)
			nodes[i] = new Node(2);
		
		nodes[0].set(1, 0.0,0.0);
		nodes[1].set(2, 1.0,0.0);
		nodes[2].set(3, 1.0,1.0);
		nodes[3].set(4, 0.0,1.0);
		nodes[4].set(5, 0.5,0.0);
		nodes[5].set(6, 1.0,0.5);
		nodes[6].set(7, 0.5,0.1);
		nodes[7].set(8, 0.0,0.5);	
		
		Element e = new Element();
		for(int i=0;i<nodes.length;i++) {
			e.addNode(nodes[i], i<4?true:false);
		}
		
		EdgeList eList = e.getEdgeList();
		for(int i=1;i<=eList.size();i++) {
			Edge edge = eList.at(i);
			System.out.println(edge);
		}		
	}
	
	public static void testRectangleNode12() {
		/*
		 * 
		 * 4--11--7--3
		 * |         |
		 * 8         10
		 * |         |
		 * 12        6
		 * |         |
		 * 1-- 5--9--2
		 * 
		 */	
		Node[] nodes = new Node[12];
		for(int i=0;i<nodes.length;i++)
			nodes[i] = new Node(2);

		double onet = 1.0/3.0;
		double twot = 2.0/3.0;
		
		nodes[0].set(  1, 0.0,0.0);
		nodes[1].set(  2, 1.0,0.0);
		nodes[2].set(  3, 1.0,1.0);
		nodes[3].set(  4, 0.0,1.0);
		nodes[4].set(  5, onet,0.0);
		nodes[5].set(  6, 1.0,onet);
		nodes[6].set(  7, twot,0.1);
		nodes[7].set(  8, 0.0,twot);	
		nodes[8].set(  9, twot,0.0);	
		nodes[9].set( 10, 1.0,twot);	
		nodes[10].set(11, onet,1.0);	
		nodes[11].set(12, 0.0,onet);	
		
		Element e = new Element();
		for(int i=0;i<nodes.length;i++) {
			e.addNode(nodes[i], i<4?true:false);
		}
		
		EdgeList eList = e.getEdgeList();
		for(int i=1;i<=eList.size();i++) {
			Edge edge = eList.at(i);
			System.out.println(edge);
		}		
	}
	
	public static void main(String[] args) {
		testTriangleNode6();
		testTriangleNode9();
		testRectangleNode8();
		testRectangleNode12();
	}
}
