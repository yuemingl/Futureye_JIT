package edu.uta.futureye.test;

import java.util.HashMap;

import org.apache.bcel.generic.ClassGen;

import edu.uta.futureye.bytecode.CompiledFunc;
import edu.uta.futureye.function.basic.FX;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.util.BytecodeUtils;
import edu.uta.futureye.util.FuncClassLoader;

public class OperatorOverloadingTest {

	public static void test1() {
		MathFunc x = FX.x;
		MathFunc y = FX.y;
		MathFunc add = x - y + 1;
		
		FuncClassLoader<CompiledFunc> fcl = new FuncClassLoader<CompiledFunc>();
		ClassGen genClass = BytecodeUtils.genClass(add, null, "add", true, false);
		CompiledFunc fadd = fcl.newInstance(genClass);
		System.out.println(fadd.apply(1.0, 2.0));
	}

	public static void test2() {
		MathFunc x = FX.x;
		MathFunc y = FX.y;
		MathFunc add = x * y;
		System.out.println(add);
		HashMap<String, MathFunc> map = new HashMap<String, MathFunc>();
		map.put(x.getVarNames().get(0), FX.r + FX.s);
		map.put(y.getVarNames().get(0), FX.r - FX.s);
		MathFunc add2 = add.compose(map);
		System.out.println(add2);
		
		FuncClassLoader<CompiledFunc> fcl = new FuncClassLoader<CompiledFunc>();
		ClassGen genClass = BytecodeUtils.genClass(add2, null, "add2", true, false);
		CompiledFunc fadd2 = fcl.newInstance(genClass);
		System.out.println(fadd2.apply(4.0, 2.0));
	}
	
	public static void test3() {
		MathFunc a = new FX("a");
		MathFunc b = new FX("b");
		MathFunc c = new FX("c");
		MathFunc x = new FX("x");
		MathFunc y = new FX("y");
		MathFunc z = new FX("z");
		MathFunc f = a*x + b*y + c*z;
		
		System.out.println(f);
		System.out.println(f.compile().apply(1,2,3,4,5,6));
		System.out.println(f.compile(new String[]{"a","b","c","x","y","z"}).apply(1,2,3,4,5,6));
	}
	
	public static void main(String[] args) {
		test1();
		test2();
		test3();
	}
}
