package edu.uta.futureye.test;

import edu.uta.futureye.core.CoordinateTransform;
import edu.uta.futureye.core.DOF;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.function.DerivativeIndicator;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.function.intf.FunctionDerivable;
import edu.uta.futureye.function.shape.SFBilinearLocal2D;
import edu.uta.futureye.function.shape.SFLinearLocal2D;

public class ShapeFunctionTest {
	

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		testSFLinearLocal2D();
		testSFBilinearLocal2D();

	}
	
	public static void testSFLinearLocal2D() {
		Node n1 = new Node(2);
		n1.set(1, 0.0,0.0);
		Node n2 = new Node(2);
		n2.set(2, 0.0,0.2);
		Node n3 = new Node(2);
		n3.set(3, 0.2,0.0);
		
		Element e = new Element();
		e.addNode(n1, true);
		e.addNode(n2, true);
		e.addNode(n3, true);
		
		SFLinearLocal2D[] shapeFun = new SFLinearLocal2D[3];
		shapeFun[0] = new SFLinearLocal2D(1);
		shapeFun[1] = new SFLinearLocal2D(2);
		shapeFun[2] = new SFLinearLocal2D(3);
		
		//Test the derivatives of shape function
		DerivativeIndicator di_x = new DerivativeIndicator("x1");
		DerivativeIndicator di_y = new DerivativeIndicator("y1");
		
		shapeFun[0].asignElement(e);
		Function SF0dx = shapeFun[0].derivative(di_x);
		Function SF0dy = shapeFun[0].derivative(di_y);
		System.out.println(SF0dx.value(null));
		System.out.println(SF0dy.value(null));
		
		shapeFun[1].asignElement(e);
		Function SF1dx = shapeFun[1].derivative(di_x);
		Function SF1dy = shapeFun[1].derivative(di_y);
		System.out.println(SF1dx.value(null));
		System.out.println(SF1dy.value(null));
		
		shapeFun[2].asignElement(e);
		Function SF2dx = shapeFun[2].derivative(di_x);
		Function SF2dy = shapeFun[2].derivative(di_y);
		System.out.println(SF2dx.value(null));
		System.out.println(SF2dy.value(null));		
		
	}
	
	public static void testSFBilinearLocal2D() {
		Node n1 = new Node(2);
		n1.set(1, -1.0,-1.0);
		Node n2 = new Node(2);
		n2.set(2, 1.0,-1.0);
		Node n3 = new Node(2);
		n3.set(3, 1.0,1.0);
		Node n4 = new Node(2);
		n4.set(4, -1.0,1.0);

		
		Element e = new Element();
		e.addNode(n1, true);
		e.addNode(n2, true);
		e.addNode(n3, true);
		e.addNode(n4, true);

		
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
			e.addDOF(j,dof);
		}
		
		//Coordinate transform and Jacbian on element e
		CoordinateTransform trans = new CoordinateTransform(2);
		trans.transformLinear2D(e);

		Variable v = new Variable();
		v.set("r", 0);
		v.set("s", 0);

		
		FunctionDerivable jac = trans.getJacobian2D();
		System.out.println(jac);
		System.out.println("jac="+jac.value(v));
		
		//Test the derivatives of shape function
		DerivativeIndicator di_x = new DerivativeIndicator("x1");
		DerivativeIndicator di_y = new DerivativeIndicator("y1");
		
		shapeFun[0].asignElement(e);
		Function SF0dx = shapeFun[0].derivative(di_x);
		Function SF0dy = shapeFun[0].derivative(di_y);
		System.out.println(SF0dx);
		System.out.println(SF0dy);
		System.out.println(SF0dx.value(v));
		
		shapeFun[1].asignElement(e);
		Function SF1dx = shapeFun[1].derivative(di_x);
		Function SF1dy = shapeFun[1].derivative(di_y);
		System.out.println(SF1dx);
		System.out.println(SF1dy);
		
		shapeFun[2].asignElement(e);
		Function SF2dx = shapeFun[2].derivative(di_x);
		Function SF2dy = shapeFun[2].derivative(di_y);
		System.out.println(SF2dx);
		System.out.println(SF2dy);	
		
		shapeFun[3].asignElement(e);
		Function SF3dx = shapeFun[3].derivative(di_x);
		Function SF3dy = shapeFun[3].derivative(di_y);
		System.out.println(SF3dx);
		System.out.println(SF3dy);		
	}

}
