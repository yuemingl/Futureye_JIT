package edu.uta.futureye.test;

import static edu.uta.futureye.function.basic.FX.x;
import static edu.uta.futureye.function.basic.FX.y;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.util.Utils;
import edu.uta.futureye.util.container.NodeList;

public class TestInterpolateOnElement {

	public static void main(String[] args) {
		System.out.println("TestInterpolateOnElement");
		NodeList nodes = new NodeList();
		nodes.add(new Node(1, -2.0,-2.0));
		nodes.add(new Node(2, 2.0,-2.0));
		nodes.add(new Node(3, 2.0,2.0));
		nodes.add(new Node(4, -2.0,2.0));
		MathFunc f = x+y*y+36;
		Element e = new Element(nodes);
		MathFunc pf = Utils.interpolateOnElement(f, e);
		System.out.println(pf.apply(0,0));
	}

}
