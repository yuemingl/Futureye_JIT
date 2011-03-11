package edu.uta.futureye.test;

import java.util.ArrayList;
import java.util.List;

import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.FAx;
import edu.uta.futureye.function.basic.FAxpb;
import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.basic.FLinear1D;
import edu.uta.futureye.function.basic.FPolynomial1D;
import edu.uta.futureye.function.basic.FX;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.function.operator.FOBasic;
import edu.uta.futureye.function.operator.FMath;

public class FunctionTest {
	
	public static void constantTest() {
		Function c0 = FC.c0;
		Function c1 = FC.c1;
		Function fx = FX.fx;
		
		System.out.println(c0.P(c0));
		System.out.println(c1.P(c1));
		System.out.println(c1.P(c0));
		System.out.println(c0.P(c1));
		System.out.println();
		System.out.println(fx.P(c0));
		System.out.println(fx.P(c1));
		System.out.println(c0.P(fx));
		System.out.println(c1.P(fx));
		System.out.println(fx.P(fx));
		System.out.println();
		
		System.out.println(c0.M(c0));
		System.out.println(c1.M(c1));
		System.out.println(c1.M(c0));
		System.out.println(c0.M(c1));
		System.out.println();
		System.out.println(fx.M(c0));
		System.out.println(fx.M(c1));
		System.out.println(c0.M(fx));
		System.out.println(c1.M(fx));
		System.out.println(fx.M(fx));
		System.out.println();
		
		System.out.println(c0.X(c0));
		System.out.println(c1.X(c1));
		System.out.println(c1.X(c0));
		System.out.println(c0.X(c1));
		System.out.println();
		System.out.println(fx.X(c0));
		System.out.println(fx.X(c1));
		System.out.println(c0.X(fx));
		System.out.println(c1.X(fx));
		System.out.println(fx.X(fx));
		System.out.println();
		
		System.out.println(c0.D(c0));
		System.out.println(c1.D(c1));
		System.out.println(c1.D(c0));
		System.out.println(c0.D(c1));
		System.out.println();
		System.out.println(fx.D(c0));
		System.out.println(fx.D(c1));
		System.out.println(c0.D(fx));
		System.out.println(c1.D(fx));
		System.out.println(fx.D(fx));
		System.out.println();
		
	}
	
	public static void test() {
		
		FAx fax1 = new FAx(3.0);
		System.out.println("fax1 = "+fax1);
		FAx fax2 = new FAx(1.0);
		System.out.println("fax2 = "+fax2);
		FAx fax3 = new FAx(0.0);
		System.out.println("fax3 = "+fax3);
		
		FAxpb faxpb1 = new FAxpb(3.0,2.0);
		System.out.println("faxpb1 = "+faxpb1);
		
		FC c0 = new FC(0.0);
		FC c1 = new FC(1.0);
		FC c2 = new FC(2.0);
		FC c3 = new FC(3.0);
		
		FX fi = FX.fx;
		Function ftest1 = FOBasic.Plus(FOBasic.Mult(c3, fi),c2);
		System.out.println("ftest1 = "+ftest1);
		Function ftest2 = FOBasic.Plus(FOBasic.Mult(c1, fi),c0);
		System.out.println("ftest2 = "+ftest2);
		Function ftest3 = FOBasic.Plus(FOBasic.Divi(c1, fi),c3);
		System.out.println("ftest3 = "+ftest3);
		
		//f(x) = c1*l1(x)+c2
		FLinear1D l1 = new FLinear1D(1,5,2,20);
		System.out.println(l1);

		Function l2 = FOBasic.Plus(FOBasic.Mult(c1, l1),c2);
		System.out.println(l2);
		Variable var = new Variable(0.5);
		System.out.println(l2.value(var));
		
		//导数：fd_x = fd'(x) , fd_xx = fd_x'(x)
//		FunctionDerivable fd = FOBasicDerivable.Plus(FOBasicDerivable.Mult(c1, fi),c3);
		Function fd = FMath.Plus(FMath.Divi(c1, fi),c3);
		Function fd_x = fd.d("x");
		Function fd_xx = fd_x.d("x");
		System.out.println("fd = "+fd);
		System.out.println("fd_x = "+fd_x);
		System.out.println("fd_xx = "+fd_xx);
		System.out.println(fd_x.value(var));
		
		//f(x) = 6*x^3 + 5*x^2 + 4*x + 3 多项式的导数
		List<Double> coef = new ArrayList<Double>();
		coef.add(3.0);
		coef.add(4.0);
		coef.add(5.0);
		coef.add(6.0);
		FPolynomial1D fp = new FPolynomial1D(coef);
		System.out.println(fp.value(new Variable(2.0)));
		Function fp_x2 = fp.d("x").d("x");
		System.out.println(fp_x2.value(new Variable(2.0)));
		Function fp_x2d = (Function)fp_x2;
		Function fp_x3 = fp_x2d.d("X");
		System.out.println(fp_x3.value(new Variable(3.0)));
		
		Function power = FOBasic.Power(c2, c3);
		System.out.println(power);
		System.out.println(power.value(null));
		
	}
	
	public static void main(String[] args) {
		
		constantTest();
		
	}

		
}
