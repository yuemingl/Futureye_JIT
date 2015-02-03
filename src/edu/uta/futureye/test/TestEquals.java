package edu.uta.futureye.test;

import edu.uta.futureye.core.Node;

public class TestEquals {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		Node node1 = new Node(100,1.0,2.0);
		Node node2 = new Node(100,1.0,2.0);
		Node node3 = new Node(101,1.5,2.0);
		Node node4 = new Node(0,1.0,2.0);
		Node node5 = new Node(0,1.6,2.0);
	
		System.out.println(node1.hashCode());
		System.out.println(node2.hashCode());
		System.out.println(node3.hashCode());
		System.out.println(node4.hashCode());
		System.out.println(node5.hashCode());

		System.out.println(node1.equals(node1));
		System.out.println(node1.equals(node2));
		System.out.println(node1.equals(node3));
		System.out.println(node1.equals(node4));
		System.out.println(node4.equals(node5));//Excepion
	}

}
