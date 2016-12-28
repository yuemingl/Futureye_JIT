package edu.uta.futureye.test;

import java.lang.reflect.Field;

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

	public static void testCompileToStaticField() {
		x.compileToStaticField(true);
		CompiledFunc cf = x.compileWithASM(new String[]{"x"});
		Field[] fds = cf.getClass().getDeclaredFields();
		System.out.println(cf.apply(10));
		for(Field f : fds) {
			try {
				System.out.println(f.getName()+"="+f.getDouble(null));
			} catch (IllegalArgumentException e) {
				e.printStackTrace();
			} catch (IllegalAccessException e) {
				e.printStackTrace();
			}
		}
		
		CompiledFunc cf2 = (x+1).compileWithASM(new String[]{"x"});
		System.out.println(cf2.apply(10));
		CompiledFunc cf3 = (x+5).compileWithASM(new String[]{"x"});
		System.out.println(cf3.apply(10));

	}
	
	public static void main(String[] args) {
//		test1();
//		test2();
		testCompileToStaticField();
	}

}
