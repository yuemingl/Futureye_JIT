package edu.uta.futureye.test;

import java.util.HashMap;
import java.util.Map;

import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.basic.FX;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.function.operator.FMath;
import edu.uta.futureye.util.container.ObjList;

public class FMathTest {
	public static void main(String[] args) {
		//u(x,y,z) = x^2+y^2+z^2
		Function u = FMath.sum(
						FX.fx.M(FX.fx),
						FX.fy.M(FX.fy),
						FX.fz.M(FX.fz));
				
		//v(x,y,x) = x*y*z
		Function v = FX.fx.M(FX.fy).M(FX.fz);
		
		//Grad(u) \cdot Grad(v)
		System.out.println(
				FMath.sum(
					u._d("x").M(v._d("x")),
					u._d("y").M(v._d("y")),
					u._d("z").M(v._d("z"))));
		
		//Grad(u) \cdot Grad(v)
		System.out.println(
				FMath.grad(u).dot(FMath.grad(v)));
		
		//invR = 1/sqrt(x^2+y^2+z^2)
		Function invR = FC.c1.D(FMath.sqrt(u));
        //Div(Grad(invR))=Laplacian(invR)
		System.out.println(FMath.div(FMath.grad(invR)));
		
//		(x + x) * y * z + (y + y) * x * z + (z + z) * x * y
//		(x + x) * y * z + (y + y) * x * z + (z + z) * x * y
//		(( - (0.5 * -0.5 * (x * x + y * y + z * z)^-1.5 * (x + x) * (x + x) + 0.5 * (x * x + y * y + z * z)^-0.5 * 2.0)) * sqrt(x * x + y * y + z * z) * sqrt(x * x + y * y + z * z) - ( - 0.5 * (x * x + y * y + z * z)^-0.5 * (x + x)) * (0.5 * (x * x + y * y + z * z)^-0.5 * (x + x) * sqrt(x * x + y * y + z * z) + sqrt(x * x + y * y + z * z) * 0.5 * (x * x + y * y + z * z)^-0.5 * (x + x))) / (sqrt(x * x + y * y + z * z) * sqrt(x * x + y * y + z * z) * sqrt(x * x + y * y + z * z) * sqrt(x * x + y * y + z * z)) + (( - (0.5 * -0.5 * (x * x + y * y + z * z)^-1.5 * (y + y) * (y + y) + 0.5 * (x * x + y * y + z * z)^-0.5 * 2.0)) * sqrt(x * x + y * y + z * z) * sqrt(x * x + y * y + z * z) - ( - 0.5 * (x * x + y * y + z * z)^-0.5 * (y + y)) * (0.5 * (x * x + y * y + z * z)^-0.5 * (y + y) * sqrt(x * x + y * y + z * z) + sqrt(x * x + y * y + z * z) * 0.5 * (x * x + y * y + z * z)^-0.5 * (y + y))) / (sqrt(x * x + y * y + z * z) * sqrt(x * x + y * y + z * z) * sqrt(x * x + y * y + z * z) * sqrt(x * x + y * y + z * z)) + (( - (0.5 * -0.5 * (x * x + y * y + z * z)^-1.5 * (z + z) * (z + z) + 0.5 * (x * x + y * y + z * z)^-0.5 * 2.0)) * sqrt(x * x + y * y + z * z) * sqrt(x * x + y * y + z * z) - ( - 0.5 * (x * x + y * y + z * z)^-0.5 * (z + z)) * (0.5 * (x * x + y * y + z * z)^-0.5 * (z + z) * sqrt(x * x + y * y + z * z) + sqrt(x * x + y * y + z * z) * 0.5 * (x * x + y * y + z * z)^-0.5 * (z + z))) / (sqrt(x * x + y * y + z * z) * sqrt(x * x + y * y + z * z) * sqrt(x * x + y * y + z * z) * sqrt(x * x + y * y + z * z))

		//invR = 1/r(x,y,z) = 1/sqrt(x^2+y^2+z^2)
		Function invR2 = FC.c1.D(FX.fr);
		Map<String, Function> fInners = 
			new HashMap<String, Function>();
		fInners.put("r", FMath.sqrt(u));
		invR2 = invR2.compose(fInners);
		
		//Grad(invR2), gradient respect to r
		System.out.println(FMath.grad(invR2));
		ObjList<String> coordSys = 
			new ObjList<String>("x","y","z");
		
        //Div_{x,y,z}(Grad_{x,y,z}(invR2))=Laplacian_{x,y,z}(invR2)
        //gradient with respect to x,y,z
        //divergence with respect to x,y,z
		System.out.println(
				FMath.div(FMath.grad(invR2, coordSys),coordSys)
				);
		
//		[-1.0 / (r * r)]
//		 ( - -1.0 * (r(x, y, z) + r(x, y, z))) / (r(x, y, z) * r(x, y, z) * r(x, y, z) * r(x, y, z)) * 0.5 * (x * x + y * y + z * z)^-0.5 * (x + x) * 0.5 * (x * x + y * y + z * z)^-0.5 * (x + x) + -1.0 / (r(x, y, z) * r(x, y, z)) * (0.5 * -0.5 * (x * x + y * y + z * z)^-1.5 * (x + x) * (x + x) + 0.5 * (x * x + y * y + z * z)^-0.5 * 2.0) + ( - -1.0 * (r(x, y, z) + r(x, y, z))) / (r(x, y, z) * r(x, y, z) * r(x, y, z) * r(x, y, z)) * 0.5 * (x * x + y * y + z * z)^-0.5 * (y + y) * 0.5 * (x * x + y * y + z * z)^-0.5 * (y + y) + -1.0 / (r(x, y, z) * r(x, y, z)) * (0.5 * -0.5 * (x * x + y * y + z * z)^-1.5 * (y + y) * (y + y) + 0.5 * (x * x + y * y + z * z)^-0.5 * 2.0) + ( - -1.0 * (r(x, y, z) + r(x, y, z))) / (r(x, y, z) * r(x, y, z) * r(x, y, z) * r(x, y, z)) * 0.5 * (x * x + y * y + z * z)^-0.5 * (z + z) * 0.5 * (x * x + y * y + z * z)^-0.5 * (z + z) + -1.0 / (r(x, y, z) * r(x, y, z)) * (0.5 * -0.5 * (x * x + y * y + z * z)^-1.5 * (z + z) * (z + z) + 0.5 * (x * x + y * y + z * z)^-0.5 * 2.0)

		
	}

}
