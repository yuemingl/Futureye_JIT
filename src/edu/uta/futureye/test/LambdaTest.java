package edu.uta.futureye.test;
import static edu.uta.futureye.function.operator.FMath.*;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.intf.MathFunc;

public class LambdaTest {

	public interface Checker {
		public double check(double a, double b);
	}
	
//	public static void main(String[] args) {
//		MathFun f = C1.A(100);
//		Checker ck = (double a, double b) -> { 
//				return f.A(a+b).apply(new Variable(50.0)); 
//			};
//		double r = ck.check(10, 2);
//		System.out.println(r);
//	}

}
