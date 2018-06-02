package edu.uta.futureye.test;

import java.lang.reflect.Field;

import edu.uta.futureye.bytecode.CompiledFunc;
import edu.uta.futureye.function.FMath;
import edu.uta.futureye.function.basic.FX;
import edu.uta.futureye.function.intf.MathFunc;
import static edu.uta.futureye.function.basic.FX.x;
import static edu.uta.futureye.function.basic.FX.y;
import static edu.uta.futureye.function.FMath.*;


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
		MathFunc fun1 = x+1;
		//set the flag to true so the result of 'fun1' will be stored to a static filed
		//after the evaluation
		fun1.compileToStaticField(true);
		CompiledFunc cf = fun1.compileWithASM(new String[]{"x"});
		System.out.println(cf.apply(10));
		//print the static field using reflection
		Field[] fds = cf.getClass().getDeclaredFields();
		for(Field f : fds) {
			try {
				System.out.println(f.getName()+"="+f.getDouble(null));
			} catch (IllegalArgumentException e) {
				e.printStackTrace();
			} catch (IllegalAccessException e) {
				e.printStackTrace();
			}
		}
		
		//Define two functions which involve 'fun1'
		//the expression of 'fun1' will not be compiled into the function 
		//instead the static field will be used as the value of 'fun1'
		CompiledFunc cf2 = (fun1+100).compileWithASM(new String[]{"x"});
		System.out.println(cf2.apply(10));
		CompiledFunc cf3 = (fun1+500).compileWithASM(new String[]{"x"});
		System.out.println(cf3.apply(10));

	}
	
	public static void testCompileToStaticField2() {
		MathFunc fun1 = x*x + y*y;
		//set the flag to true so the result of 'fun1' will be stored to a static filed
		//after the evaluation
		fun1.compileToStaticField(true);
		CompiledFunc cf = fun1.compileWithASM(new String[]{"x","y"});
		System.out.println(cf.apply(3, 4));
		//print the static field using reflection
		Field[] fds = cf.getClass().getDeclaredFields();
		for(Field f : fds) {
			try {
				System.out.println(f.getName()+"="+f.getDouble(null));
			} catch (IllegalArgumentException e) {
				e.printStackTrace();
			} catch (IllegalAccessException e) {
				e.printStackTrace();
			}
		}
		
		//Define two functions which involve 'fun1'
		//the expression of 'fun1' will not be compiled into the function 
		//instead the static field will be used as the value of 'fun1'
		CompiledFunc cf2 = (sqrt(fun1)).compileWithASM(new String[]{"x"});
		System.out.println(cf2.apply(10));
		CompiledFunc cf3 = (sqrt(fun1)).compileWithASM(new String[]{"x"});
		System.out.println(cf3.apply(10));

	}
	
	
	public static void main(String[] args) {
		test1();
//		test2();
//		testCompileToStaticField();
		testCompileToStaticField2();
	}

}
