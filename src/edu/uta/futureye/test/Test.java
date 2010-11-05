package edu.uta.futureye.test;

import java.util.LinkedList;
import java.util.List;

import edu.uta.futureye.util.Utils;

public class Test {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
//		Variable v = new Variable();
//		v.add("x", 3.0);
//		System.out.println(v.getValue("y"));
		
		List<String> a = new LinkedList<String>();
		a.add("x");
		a.add("z");
		List<String> b = new LinkedList<String>();
		b.add("x");
		b.add("y");
		List<String> c= Utils.mergeList(a, b);
		System.out.println(c);
		
//		double eps = 1e-6;
//		double pow = Math.pow(Math.E, -0.25*0.0001/eps);
//		double delta = 0.5*pow/Math.sqrt(Math.PI*eps);
//		System.out.println(delta);
		
		System.out.println("a".compareTo(null));
	}

}
