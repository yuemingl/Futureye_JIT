package edu.uta.futureye.test;

import java.util.HashMap;
import java.util.Map;
import edu.uta.futureye.algebra.SpaceVector;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.util.container.ObjList;
import static edu.uta.futureye.function.operator.FMath.*;

public class FMathTest {
	
	public static void test1() {
		//u(x,y,z) = x^2+y^2+z^2
		MathFunc u = sum(X.M(X), Y.M(Y), Z.M(Z));

		//v(x,y,x) = x*y*z
		MathFunc v = X.M(Y).M(Z);
		
		//Grad(u) \cdot Grad(v)
		System.out.println(
				sum(
					u._d("x").M(v._d("x")),
					u._d("y").M(v._d("y")),
					u._d("z").M(v._d("z"))));
		
		//Grad(u) \cdot Grad(v)
		System.out.println(grad(u).dot(grad(v)));
		
		//invR = 1/sqrt(x^2+y^2+z^2)
		MathFunc invR = C1.D(sqrt(u));
		//Div(Grad(invR))=Laplacian(invR)
		System.out.println(div(grad(invR)));
//		(x + x) * y * z + (y + y) * x * z + (z + z) * x * y
//		(x + x) * y * z + (y + y) * x * z + (z + z) * x * y
//		(( - (0.5 * -0.5 * (x * x + y * y + z * z)^-1.5 * (x + x) * (x + x) + 0.5 * (x * x + y * y + z * z)^-0.5 * 2.0)) * sqrt(x * x + y * y + z * z) * sqrt(x * x + y * y + z * z) - ( - 0.5 * (x * x + y * y + z * z)^-0.5 * (x + x)) * (0.5 * (x * x + y * y + z * z)^-0.5 * (x + x) * sqrt(x * x + y * y + z * z) + sqrt(x * x + y * y + z * z) * 0.5 * (x * x + y * y + z * z)^-0.5 * (x + x))) / (sqrt(x * x + y * y + z * z) * sqrt(x * x + y * y + z * z) * sqrt(x * x + y * y + z * z) * sqrt(x * x + y * y + z * z)) + (( - (0.5 * -0.5 * (x * x + y * y + z * z)^-1.5 * (y + y) * (y + y) + 0.5 * (x * x + y * y + z * z)^-0.5 * 2.0)) * sqrt(x * x + y * y + z * z) * sqrt(x * x + y * y + z * z) - ( - 0.5 * (x * x + y * y + z * z)^-0.5 * (y + y)) * (0.5 * (x * x + y * y + z * z)^-0.5 * (y + y) * sqrt(x * x + y * y + z * z) + sqrt(x * x + y * y + z * z) * 0.5 * (x * x + y * y + z * z)^-0.5 * (y + y))) / (sqrt(x * x + y * y + z * z) * sqrt(x * x + y * y + z * z) * sqrt(x * x + y * y + z * z) * sqrt(x * x + y * y + z * z)) + (( - (0.5 * -0.5 * (x * x + y * y + z * z)^-1.5 * (z + z) * (z + z) + 0.5 * (x * x + y * y + z * z)^-0.5 * 2.0)) * sqrt(x * x + y * y + z * z) * sqrt(x * x + y * y + z * z) - ( - 0.5 * (x * x + y * y + z * z)^-0.5 * (z + z)) * (0.5 * (x * x + y * y + z * z)^-0.5 * (z + z) * sqrt(x * x + y * y + z * z) + sqrt(x * x + y * y + z * z) * 0.5 * (x * x + y * y + z * z)^-0.5 * (z + z))) / (sqrt(x * x + y * y + z * z) * sqrt(x * x + y * y + z * z) * sqrt(x * x + y * y + z * z) * sqrt(x * x + y * y + z * z))

		//invR = 1/r(x,y,z) = 1/sqrt(x^2+y^2+z^2)
		MathFunc invR2 = C1.D(R);
		Map<String, MathFunc> fInners = 
			new HashMap<String, MathFunc>();
		fInners.put("r", sqrt(u));
		invR2 = invR2.compose(fInners);
		
		//Grad(invR2), gradient respect to r
		System.out.println(grad(invR2));
		ObjList<String> coordSys = 
			new ObjList<String>("x","y","z");
		
		//Div_{x,y,z}(Grad_{x,y,z}(invR2))=Laplacian_{x,y,z}(invR2)
		//gradient with respect to x,y,z
		//divergence with respect to x,y,z
		System.out.println(div(grad(invR2, coordSys),coordSys));
//		[-1.0 / (r * r)]
//		 ( - -1.0 * (r(x, y, z) + r(x, y, z))) / (r(x, y, z) * r(x, y, z) * r(x, y, z) * r(x, y, z)) * 0.5 * (x * x + y * y + z * z)^-0.5 * (x + x) * 0.5 * (x * x + y * y + z * z)^-0.5 * (x + x) + -1.0 / (r(x, y, z) * r(x, y, z)) * (0.5 * -0.5 * (x * x + y * y + z * z)^-1.5 * (x + x) * (x + x) + 0.5 * (x * x + y * y + z * z)^-0.5 * 2.0) + ( - -1.0 * (r(x, y, z) + r(x, y, z))) / (r(x, y, z) * r(x, y, z) * r(x, y, z) * r(x, y, z)) * 0.5 * (x * x + y * y + z * z)^-0.5 * (y + y) * 0.5 * (x * x + y * y + z * z)^-0.5 * (y + y) + -1.0 / (r(x, y, z) * r(x, y, z)) * (0.5 * -0.5 * (x * x + y * y + z * z)^-1.5 * (y + y) * (y + y) + 0.5 * (x * x + y * y + z * z)^-0.5 * 2.0) + ( - -1.0 * (r(x, y, z) + r(x, y, z))) / (r(x, y, z) * r(x, y, z) * r(x, y, z) * r(x, y, z)) * 0.5 * (x * x + y * y + z * z)^-0.5 * (z + z) * 0.5 * (x * x + y * y + z * z)^-0.5 * (z + z) + -1.0 / (r(x, y, z) * r(x, y, z)) * (0.5 * -0.5 * (x * x + y * y + z * z)^-1.5 * (z + z) * (z + z) + 0.5 * (x * x + y * y + z * z)^-0.5 * 2.0)
	}
	
	public static void test2() {
		double[] dataSet = {2,4,5,5,7};
		SpaceVector v = new SpaceVector(dataSet);
		System.out.println(mean(v));
		System.out.println(variance(v));
		System.out.println(SD(v));
		System.out.println(averageAbsoluteDeviation(v));
	}
	
	public static void main(String[] args) {
		//test1();
		test2();
	}

}
