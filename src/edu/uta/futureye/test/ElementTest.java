package edu.uta.futureye.test;

import edu.uta.futureye.core.EdgeLocal;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.util.list.NodeList;
import edu.uta.futureye.util.list.ObjList;

public class ElementTest {
	
	public static void testTriangleNode3_Clockwise() {
		/*
		 * À≥ ±’Î≤‚ ‘
		 *  2
		 *  | \
		 *  |  \
		 *  1---3 
		 */
		NodeList nodes = new NodeList();
		nodes.add(new Node(1, 0.0,0.0));
		nodes.add(new Node(2, 0.0,0.2));
		nodes.add(new Node(3, 0.2,0.0));
		
		Element e = new Element(nodes);
		
		System.out.println("TriangleNode3: \n"+e.vertices());
		System.out.println(e.nodes);
		e.adjustVerticeToCounterClockwise();
		System.out.println("adjustVerticeToCounterClockwise: \n"+e.vertices());
		System.out.println(e.nodes);
	}

	public static void testTriangleNode3() {
		/*
		 *  3
		 *  | \
		 *  |  \
		 *  1---2 
		 */
		NodeList nodes = new NodeList();
		nodes.add(new Node(1, 0.0,0.0));
		nodes.add(new Node(2, 0.2,0.0));
		nodes.add(new Node(3, 0.0,0.2));
		
		Element e = new Element(nodes);
		
		System.out.println("TriangleNode3:");
		e.updateJacobinLinear2D();
		System.out.println(e.getJacobin());
		
		ObjList<EdgeLocal> eList = e.edges();
		for(int i=1;i<=eList.size();i++) {
			EdgeLocal edge = eList.at(i);
			System.out.println(edge);
		}
	}
	

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
		NodeList nodes = new NodeList();
		nodes.add(new Node(1, 0.0,0.0));
		nodes.add(new Node(2, 1.0,0.0));
		nodes.add(new Node(3, 0.0,1.0));
		nodes.add(new Node(4, 0.5,0.0));
		nodes.add(new Node(5, 0.5,0.5));
		nodes.add(new Node(6, 0.0,0.5));
		
		Element e = new Element(nodes);
		
		System.out.println("TriangleNode6:");
		e.updateJacobinLinear2D();
		System.out.println(e.getJacobin());
		
		ObjList<EdgeLocal> eList = e.edges();
		for(int i=1;i<=eList.size();i++) {
			EdgeLocal edge = eList.at(i);
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
		NodeList list = new NodeList();
		list.fromArray(nodes);
		System.out.println("TriangleNode9:");
		Element e = new Element(list);
		
		ObjList<EdgeLocal> eList = e.edges();
		for(int i=1;i<=eList.size();i++) {
			EdgeLocal edge = eList.at(i);
			System.out.println(edge);
		}
	}	
	
	
	public static void testRectangleNode4_Clockwise() {
		/*
		 * À≥ ±’Î≤‚ ‘
		 * 2----3
		 * |    |
		 * |    |
		 * 1----4
		 * 
		 */	
		NodeList nodes = new NodeList();
		nodes.add(new Node(1, -1.0,-1.0));
		nodes.add(new Node(2, -1.0,1.0));
		nodes.add(new Node(3, 1.0,1.0));
		nodes.add(new Node(4, 1.0,-1.0));
		
		Element e = new Element(nodes);
		
		System.out.println("RectangleNode4: \n"+e.vertices());
		System.out.println(e.nodes);
		e.adjustVerticeToCounterClockwise();
		System.out.println("adjustVerticeToCounterClockwise: \n"+e.vertices());
		System.out.println(e.nodes);
	}	
	
	public static void testRectangleNode4() {
		/*
		 * 4----3
		 * |    |
		 * |    |
		 * 1----2
		 * 
		 */	
		NodeList nodes = new NodeList();
		nodes.add(new Node(1, -1.0,-1.0));
		nodes.add(new Node(2, 1.0,-1.0));
		nodes.add(new Node(3, 1.0,1.0));
		nodes.add(new Node(4, -1.0,1.0));
		
		System.out.println("RectangleNode4:");
		Element e = new Element(nodes);
		e.updateJacobinLinear2D();
		System.out.println(e.getJacobin());
		
		ObjList<EdgeLocal> eList = e.edges();
		for(int i=1;i<=eList.size();i++) {
			EdgeLocal edge = eList.at(i);
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
		NodeList list = new NodeList();
		list.fromArray(nodes);
		System.out.println("RectangleNode8:");
		Element e = new Element(list);
		
		ObjList<EdgeLocal> eList = e.edges();
		for(int i=1;i<=eList.size();i++) {
			EdgeLocal edge = eList.at(i);
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
		NodeList list = new NodeList();
		list.fromArray(nodes);
		Element e = new Element(list);
		
		System.out.println("RectangleNode12:");
		ObjList<EdgeLocal> eList = e.edges();
		for(int i=1;i<=eList.size();i++) {
			EdgeLocal edge = eList.at(i);
			System.out.println(edge);
		}	
	}
	
	public static void main(String[] args) {
		testTriangleNode3_Clockwise();
		testRectangleNode4_Clockwise();
//		testTriangleNode3();
//		testRectangleNode4();
//		testTriangleNode6();
//		testTriangleNode9();
//		testRectangleNode8();
//		testRectangleNode12();
	}
}
