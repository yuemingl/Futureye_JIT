package edu.uta.futureye.test;

import java.util.HashMap;

import com.sun.org.apache.bcel.internal.generic.ClassGen;

import edu.uta.futureye.bytecode.CompiledFunc;
import edu.uta.futureye.function.basic.FX;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.util.BytecodeUtils;
import edu.uta.futureye.util.FuncClassLoader;

public class BytecodeTest {

	public static void test1() {
		MathFunc x = FX.fx;
		MathFunc y = FX.fy;
		MathFunc add = x.S(y);
		
		FuncClassLoader<CompiledFunc> fcl = new FuncClassLoader<CompiledFunc>();
		ClassGen genClass = BytecodeUtils.genClass(add, "add", true, false);
		CompiledFunc fadd = fcl.newInstance(genClass);
		System.out.println(fadd.apply(1.0, 2.0));
	}

	public static void test2() {
		MathFunc x = FX.fx;
		MathFunc y = FX.fy;
		MathFunc add = x.M(y);
		System.out.println(add);
		HashMap<String, MathFunc> map = new HashMap<String, MathFunc>();
		map.put(x.getVarNames().get(0), FX.fr.A(FX.fs));
		map.put(y.getVarNames().get(0), FX.fr.S(FX.fs));
		MathFunc add2 = add.compose(map);
		System.out.println(add2);
		
		FuncClassLoader<CompiledFunc> fcl = new FuncClassLoader<CompiledFunc>();
		ClassGen genClass = BytecodeUtils.genClass(add2, "add2", true, false);
		CompiledFunc fadd2 = fcl.newInstance(genClass);
		System.out.println(fadd2.apply(4.0, 2.0));
	}
	
	public static void main(String[] args) {
		test1();
		test2();
	}

}
