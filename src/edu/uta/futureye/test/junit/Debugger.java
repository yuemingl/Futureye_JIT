package edu.uta.futureye.test.junit;

import static edu.uta.futureye.function.FMath.x;
import static edu.uta.futureye.function.FMath.y;
import edu.uta.futureye.function.FMath;

public class Debugger {

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		System.out.println(FMath.pow(x,y).compileWithASM(new String[]{"x","y"}).apply(2.0, 5));

	}

}
