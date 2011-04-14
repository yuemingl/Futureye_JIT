package edu.uta.futureye.test;

import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.function.intf.ScalarShapeFunction;
import edu.uta.futureye.function.operator.FOIntegrate;
import edu.uta.futureye.lib.shapefun.SFBilinearLocal2D;
import edu.uta.futureye.lib.shapefun.SFLinearLocal1D;
import edu.uta.futureye.lib.shapefun.SFLinearLocal2D;

public class IntegrationTest {
	
	public static void testSFLinearLocal1D() {
		ScalarShapeFunction sf1d1 = new SFLinearLocal1D(1);
		ScalarShapeFunction sf1d2 = new SFLinearLocal1D(2);
		
		Function integrand = sf1d1.A(sf1d2);
		
		Function integral;

		integral = FOIntegrate.intOnTriangleRefElement(sf1d1, 4);
		System.out.println(integral.value(null));
		
		integral = FOIntegrate.intOnTriangleRefElement(sf1d2, 4);
		System.out.println(integral.value(null));
		
		integral = FOIntegrate.intOnLinearRefElement(integrand, 4);
		System.out.println(integral.value(null));
		
//		0.16666666666666669
//		0.33333333333333337
//		2.0
	}
	
	public static void testSFLinearLocal2D() {
		SFLinearLocal2D[] shapeFun = new SFLinearLocal2D[3];
		shapeFun[0] = new SFLinearLocal2D(1);
		shapeFun[1] = new SFLinearLocal2D(2);
		shapeFun[2] = new SFLinearLocal2D(3);
		
		Function integrand = shapeFun[0].A(shapeFun[1]);
		integrand = integrand.A(shapeFun[2]);
		
		Function integral;
		
		integral = FOIntegrate.intOnTriangleRefElement(shapeFun[0], 3);
		System.out.println(integral.value(null));
		
		integral = FOIntegrate.intOnTriangleRefElement(shapeFun[1], 3);
		System.out.println(integral.value(null));
		
		integral = FOIntegrate.intOnTriangleRefElement(shapeFun[2], 3);
		System.out.println(integral.value(null));
		
		integral = FOIntegrate.intOnTriangleRefElement(integrand, 3);
		System.out.println(integral.value(null));
		
//		0.16666666666666666
//		0.16666666666666666
//		0.16666666666666666
//		0.5		
	}
	
	public static void testSFBilinearLocal2D() {
		SFBilinearLocal2D[] shapeFun = new SFBilinearLocal2D[4];
		shapeFun[0] = new SFBilinearLocal2D(1);
		shapeFun[1] = new SFBilinearLocal2D(2);
		shapeFun[2] = new SFBilinearLocal2D(3);
		shapeFun[3] = new SFBilinearLocal2D(3);
		
		Function integrand = shapeFun[0].A(shapeFun[1]);
		integrand = integrand.A(shapeFun[2]);
		integrand = integrand.A(shapeFun[3]);
		
		Function integral;
		
		integral = FOIntegrate.intOnTriangleRefElement(shapeFun[0], 2);
		System.out.println(integral.value(null));
		
		integral = FOIntegrate.intOnTriangleRefElement(shapeFun[1], 2);
		System.out.println(integral.value(null));
		
		integral = FOIntegrate.intOnTriangleRefElement(shapeFun[2], 2);
		System.out.println(integral.value(null));
		
		integral = FOIntegrate.intOnTriangleRefElement(shapeFun[3], 2);
		System.out.println(integral.value(null));
		
		integral = FOIntegrate.intOnTriangleRefElement(integrand, 2);
		System.out.println(integral.value(null));
		
//		0.055555555555555566
//		0.11111111111111112
//		0.2222222222222222
//		0.2222222222222222
//		0.6111111111111112		
	}
	
	public static void main(String[] args) {
		testSFLinearLocal1D();
		testSFLinearLocal2D();
		testSFBilinearLocal2D();
		
	}
}
