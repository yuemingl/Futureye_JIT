package edu.uta.futureye.test;

import edu.uta.futureye.core.Node;
import edu.uta.futureye.util.Constant;
import edu.uta.futureye.util.Utils;

public class UtilsTest {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		double []o = {0.0,0.0,0.0};
		double []a = {1.0,0.0,0.0};
		double []b = {0.0,1.0,0.0};
		double []c = {1.0,1.0,1.0};
		
		double area = Utils.getSphereTriangleArea(1,
				new Node().set(0, o),
				new Node().set(0, a),
				new Node().set(0, b),
				new Node().set(0, c)
				);
		
		System.out.println(area);
		if(Math.abs(0.5235987755982987-area)>Constant.eps)
			System.out.println("Test Fail!!!");

	}

}
