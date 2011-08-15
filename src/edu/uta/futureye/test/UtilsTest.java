package edu.uta.futureye.test;

import edu.uta.futureye.core.Node;
import edu.uta.futureye.core.Vertex;
import edu.uta.futureye.core.geometry.Point;
import edu.uta.futureye.util.Constant;
import edu.uta.futureye.util.Utils;

public class UtilsTest {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		Point []p = new Vertex[4];
//		p[0] = new Vertex(1, 1.0, 0.5);
//		p[1] = new Vertex(2, 1.5, 1.0);
//		p[2] = new Vertex(3, 1.0, 1.5);
//		p[3] = new Vertex(4, 0.5, 1.0);
		
		p[0] = new Vertex(1, 0.0, 2.75);
		p[1] = new Vertex(2, 0.25, 2.75);
		p[2] = new Vertex(3, 0.25, 3.0);
		p[3] = new Vertex(4, 0.0, 3.0);		
		double[] f = {2.0,3.0,4.0,10.0};
		
//		p[0] = new Vertex(1, 0.0, 3.0);		
//		p[1] = new Vertex(2, 0.0, 2.75);
//		p[2] = new Vertex(3, 0.25, 2.75);
//		p[3] = new Vertex(4, 0.25, 3.0);
//		double[] f = {10.0,2.0,3.0,4.0};
		
		double[] a = Utils.computeBilinearFunctionCoef(p, f);
		
		for(int i=0;i<4;i++) {
			double x = p[i].coord(1);
			double y = p[i].coord(2);
			System.out.format("1 %f %f %f \r\n", x,y,x*y);
		}
		
		for(int i=0;i<4;i++) {
			System.out.print(a[i]+" ");
		}
		System.out.println();
		for(int i=0;i<4;i++) {
			double x = p[i].coord(1);
			double y = p[i].coord(2);
			double fv = a[0]+a[1]*x+a[2]*y+a[3]*x*y;
			System.out.println(fv);
		}

	}
	
	public static void test1() {
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
