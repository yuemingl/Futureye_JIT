package edu.uta.futureye.function;

import edu.uta.futureye.function.intf.MathFunc;

public class FTest extends AbstractMathFunc {

	@Override
	public double apply(double... args) {
		// TODO Auto-generated method stub
		return 0;
	}

	public static void main(String[]args) {
		System.out.println(new FTest());
		
		MathFunc f = new AbstractMathFunc("x","y") {
	
			@Override
			public double apply(double... args) {
				// TODO Auto-generated method stub
				return 0;
			}
	
		};
		System.out.println(f);
	}

}
