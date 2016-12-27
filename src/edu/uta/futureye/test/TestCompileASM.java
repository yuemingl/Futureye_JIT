package edu.uta.futureye.test;

import edu.uta.futureye.bytecode.CompiledFunc;
import edu.uta.futureye.function.basic.FX;
import static edu.uta.futureye.function.basic.FX.x;

public class TestCompileASM {

	public static void test1() {
		CompiledFunc cf = FX.x.compileWithASM(new String[]{"x"});
		//CompiledFunc cf = FX.x.compile(new String[]{"x"});
		System.out.println(cf.apply(10));
	}

	public static void test2() {
		CompiledFunc cf = (x+1).compileWithASM(new String[]{"x"});
		System.out.println(cf.apply(10));
	}

	
	public static void main(String[] args) {
		test1();
		test2();
	}

}
