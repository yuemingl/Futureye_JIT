package edu.uta.futureye.test;

import edu.uta.futureye.core.CoordinateTransform;
import edu.uta.futureye.core.DOF;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.lib.shapefun.SFBilinearLocal2D;
import edu.uta.futureye.lib.shapefun.SFLinearLocal2DTest;
import edu.uta.futureye.util.container.NodeList;

public class ShapeFunctionTest {
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		testSFLinearLocal2D();
		testSFBilinearLocal2D();

	}
	
	public static void testSFLinearLocal2D() {
		NodeList nodes = new NodeList();
		nodes.add(new Node(1, 0.0,0.0));
		nodes.add(new Node(2, 0.2,0.0));
		nodes.add(new Node(3, 0.0,0.2));
		
		Element e = new Element(nodes);
		
//		SFLinearLocal2D[] shapeFun = new SFLinearLocal2D[3];
//		shapeFun[0] = new SFLinearLocal2D(1);
//		shapeFun[1] = new SFLinearLocal2D(2);
//		shapeFun[2] = new SFLinearLocal2D(3);
		SFLinearLocal2DTest[] shapeFun = new SFLinearLocal2DTest[3];
		shapeFun[0] = new SFLinearLocal2DTest(1);
		shapeFun[1] = new SFLinearLocal2DTest(2);
		shapeFun[2] = new SFLinearLocal2DTest(3);
		
		//Test the derivatives of shape function
		
		shapeFun[0].asignElement(e);
		System.out.println(shapeFun[0]);
		Function SF0dx = shapeFun[0]._d("x");
		Function SF0dy = shapeFun[0]._d("y");
		System.out.println(SF0dx);
		System.out.println(SF0dx.value(null));
		System.out.println(SF0dy);
		System.out.println(SF0dy.value(null));
		
		shapeFun[1].asignElement(e);
		System.out.println(shapeFun[1]);
		Function SF1dx = shapeFun[1]._d("x");
		Function SF1dy = shapeFun[1]._d("y");
		System.out.println(SF1dx.value(null));
		System.out.println(SF1dy.value(null));
		
		shapeFun[2].asignElement(e);
		System.out.println(shapeFun[2]);
		Function SF2dx = shapeFun[2]._d("x");
		Function SF2dy = shapeFun[2]._d("y");
		System.out.println(SF2dx.value(null));
		System.out.println(SF2dy.value(null));		
		
	}
	
	public static void testSFBilinearLocal2D() {
		NodeList nodes = new NodeList();
		nodes.add(new Node(1, -1.0,-1.0));
		nodes.add(new Node(2, 1.0,-1.0));
		nodes.add(new Node(3, 1.0,1.0));
		nodes.add(new Node(4, -1.0,1.0));
		
		Element e = new Element(nodes);

		SFBilinearLocal2D[] shapeFun = new SFBilinearLocal2D[4];
		shapeFun[0] = new SFBilinearLocal2D(1);
		shapeFun[1] = new SFBilinearLocal2D(2);
		shapeFun[2] = new SFBilinearLocal2D(3);
		shapeFun[3] = new SFBilinearLocal2D(4);
		Variable v0 = new Variable();
		v0.set("r", 1.0);
		v0.set("s", 0.0);
		System.out.println(shapeFun[0].value(v0));
		
		//Asign degree of freedom to nodes
		for(int j=1;j<=e.nodes.size();j++) {
			//Asign shape function to DOF
			DOF dof = new DOF(j,e.nodes.at(j).globalIndex,shapeFun[j-1]);
			e.addNodeDOF(j,dof);
		}
		
		//Coordinate transform and Jacbian on element e
		CoordinateTransform trans = new CoordinateTransform(2);
		trans.transformLinear2D(e);

		Variable v = new Variable();
		v.set("r", 0);
		v.set("s", 0);

		
		Function jac = trans.getJacobian2D();
		System.out.println(jac);
		System.out.println("jac="+jac.value(v));
		
		//Test the derivatives of shape function
		
		shapeFun[0].asignElement(e);
		Function SF0dx = shapeFun[0]._d("x");
		Function SF0dy = shapeFun[0]._d("y");
		System.out.println(SF0dx);
		System.out.println("SF0dx("+v+")="+SF0dx.value(v));
		System.out.println(SF0dy);
		System.out.println("SF0dy("+v+")="+SF0dy.value(v));
		
		shapeFun[1].asignElement(e);
		Function SF1dx = shapeFun[1]._d("x");
		Function SF1dy = shapeFun[1]._d("y");
		System.out.println(SF1dx);
		System.out.println(SF1dy);
		
		shapeFun[2].asignElement(e);
		Function SF2dx = shapeFun[2]._d("x");
		Function SF2dy = shapeFun[2]._d("y");
		System.out.println(SF2dx);
		System.out.println(SF2dy);	
		
		shapeFun[3].asignElement(e);
		Function SF3dx = shapeFun[3]._d("x");
		Function SF3dy = shapeFun[3]._d("y");
		System.out.println(SF3dx);
		System.out.println(SF3dy);		
	}

}
