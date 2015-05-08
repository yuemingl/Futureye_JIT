package edu.uta.futureye.test;

import java.util.HashMap;

import com.sun.org.apache.bcel.internal.generic.ClassGen;

import edu.uta.futureye.bytecode.CompiledFunc;
import edu.uta.futureye.function.FSqrt;
import edu.uta.futureye.function.basic.FAx;
import edu.uta.futureye.function.basic.FAxpb;
import edu.uta.futureye.function.basic.FCos;
import edu.uta.futureye.function.basic.FSin;
import edu.uta.futureye.function.basic.FTan;
import edu.uta.futureye.function.basic.FX;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.function.operator.FMath;
import edu.uta.futureye.util.BytecodeUtils;
import edu.uta.futureye.util.FuncClassLoader;

public class BytecodeTest {

	public static void test1() {
		MathFunc x = FX.x;
		MathFunc y = FX.y;
		MathFunc add = x.S(y);
		
		FuncClassLoader<CompiledFunc> fcl = new FuncClassLoader<CompiledFunc>();
		ClassGen genClass = BytecodeUtils.genClass(add, "add", true, false);
		CompiledFunc fadd = fcl.newInstance(genClass);
		System.out.println(fadd.apply(1.0, 2.0));
	}

	public static void test2() {
		MathFunc x = FX.x;
		MathFunc y = FX.y;
		MathFunc add = x.M(y);
		System.out.println(add);
		HashMap<String, MathFunc> map = new HashMap<String, MathFunc>();
		map.put(x.getVarNames().get(0), FX.r.A(FX.s));
		map.put(y.getVarNames().get(0), FX.r.S(FX.s));
		MathFunc add2 = add.compose(map);
		System.out.println(add2);
		
		FuncClassLoader<CompiledFunc> fcl = new FuncClassLoader<CompiledFunc>();
		ClassGen genClass = BytecodeUtils.genClass(add2, "add2", true, false);
		CompiledFunc fadd2 = fcl.newInstance(genClass);
		System.out.println(fadd2.apply(4.0, 2.0));
	}
	
	public static void test3() {
		FAx fax = new FAx(2.0);
		fax.setName("fax");
		CompiledFunc cfax = fax.compile();
		System.out.println(cfax.apply(0.1));
	}

	public static void test4() {
		FSin sin = new FSin("x");
		FCos cos = new FCos("y");
		MathFunc f = sin * cos;
		System.out.println(f);
		System.out.println(f.getVarNames());
		System.out.println(f.apply(Math.PI/2, Math.PI/4));
		
		CompiledFunc cf = f.compile();
		System.out.println(cf.apply(Math.PI/2, Math.PI/4));
	}

	public static void test5() {
		MathFunc f = FX.x * FX.y * FX.z + 1;
		System.out.println(f);
		System.out.println(f.getVarNames());
		System.out.println(f.apply(2,3,4));
		CompiledFunc cf = f.compile();
		System.out.println(cf.apply(2,3,4));
	}
	
	public static void test6() {
		FTan tan = new FTan("x");
		System.out.println(tan.diff("x"));
		
		FAx fax = new FAx("x",2.0);
		System.out.println(fax);
		System.out.println(fax.diff("x"));
		System.out.println(fax.compile().apply(5.0));
		
		FAxpb faxpb = new FAxpb("x",2.0, 3.0);
		System.out.println(faxpb);
		System.out.println(faxpb.diff("x"));
		System.out.println(faxpb.compile().apply(5.0));
	}
	
	public static void test7() {
		FSin sin = new FSin();
		MathFunc fun = FMath.pow(sin, 2);
		System.out.println(fun.toString());
		System.out.println(fun.diff("x"));
		
		System.out.println(FMath.sqrt(sin));
		System.out.println(FMath.sqrt(sin).diff("x"));
		
		FSqrt sqrt = new FSqrt();
		System.out.println(sqrt.compile().apply(16));
		System.out.println(FMath.sqrt(FX.x).compile().apply(16));
		
	}
	public static void main(String[] args) {
		test1();
		test2();
		test3();
		test4();
		test5();
		test6();
		test7();
	}

}
