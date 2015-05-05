package edu.uta.futureye.test;

import java.util.HashMap;

import com.sun.org.apache.bcel.internal.generic.ClassGen;

import edu.uta.futureye.bytecode.BytecodeFunc;
import edu.uta.futureye.function.basic.FX;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.util.BytecodeUtils;
import edu.uta.futureye.util.FuncClassLoader;

public class OperatorOverloadingTest {

	public static void test1() {
		MathFunc x = FX.fx;
		MathFunc y = FX.fy;
		MathFunc add = x - y;
		
		FuncClassLoader<BytecodeFunc> fcl = new FuncClassLoader<BytecodeFunc>();
		ClassGen genClass = BytecodeUtils.genClassBytecodeFunc(add, "add", true, false);
		BytecodeFunc fadd = fcl.newInstance(genClass);
		System.out.println(fadd.apply(1.0, 2.0));
	}

	public static void test2() {
		MathFunc x = FX.fx;
		MathFunc y = FX.fy;
		MathFunc add = x * y;
		System.out.println(add);
		HashMap<String, MathFunc> map = new HashMap<String, MathFunc>();
		map.put(x.getVarNames().get(0), FX.fr + FX.fs);
		map.put(y.getVarNames().get(0), FX.fr - FX.fs);
		MathFunc add2 = add.compose(map);
		System.out.println(add2);
		
		FuncClassLoader<BytecodeFunc> fcl = new FuncClassLoader<BytecodeFunc>();
		ClassGen genClass = BytecodeUtils.genClassBytecodeFunc(add2, "add2", true, false);
		BytecodeFunc fadd2 = fcl.newInstance(genClass);
		System.out.println(fadd2.apply(4.0, 2.0));
	}
	
	public static void main(String[] args) {
		test1();
		test2();
	}
}
