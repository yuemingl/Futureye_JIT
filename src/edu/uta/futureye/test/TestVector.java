package edu.uta.futureye.test;

import edu.uta.futureye.algebra.Vector;

public class TestVector {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		Vector v = new Vector(3);
		v.set(1, 1.0);
		v.set(2, 2.0);
		v.set(3, 3.0);
		v = Vector.axmy(1.0, v, v);
		v.print();
	}

}
