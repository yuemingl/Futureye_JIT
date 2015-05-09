package edu.uta.futureye.function;

import edu.uta.futureye.function.intf.MathFunc;

public class FTest extends AbstractMathFunc {


public static void main(String[]args) {
	System.out.println(new FTest());
	
	MathFunc f = new AbstractMathFunc("x","y") {

		@Override
		public double apply(double... args) {
			// TODO Auto-generated method stub
			return 0;
		}

		@Override
		public double apply(Variable v) {
			// TODO Auto-generated method stub
			return 0;
		}
		
	};
	System.out.println(f);
}

@Override
public double apply(double... args) {
	// TODO Auto-generated method stub
	return 0;
}

@Override
public double apply(Variable v) {
	// TODO Auto-generated method stub
	return 0;
}

}
