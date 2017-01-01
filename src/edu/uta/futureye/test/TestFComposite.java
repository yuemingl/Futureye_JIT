package edu.uta.futureye.test;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import edu.uta.futureye.bytecode.CompiledFunc;
import edu.uta.futureye.function.basic.FX;
import edu.uta.futureye.function.intf.MathFunc;
import static edu.uta.futureye.function.FMath.*;

public class TestFComposite {
	public static void check(String a, String b) {
		if(a.equals(b))
			System.out.println("pass");
		else 
			System.out.println("!!!FAIL!!!   " + a + " != " + b);
	}
	
	public static void test0() {
		MathFunc f = r*s + r + s + 1;
		Map<String, MathFunc> fInners = new HashMap<String, MathFunc>();
		fInners.put("r", x*x);
		fInners.put("s", y+1);
		MathFunc fc = f.compose(fInners);
		
		String fc_expr = fc.toString();
		System.out.println(fc.compile().apply(new double[]{2.0,3.0})); //25.0
		System.out.println(fc.compile(new String[]{"a","b","x","y"}).apply(new double[]{1.0,2.0,2.0,3.0}));
		check(fc_expr, fc.toString()); //f(x,y) = (x*x)*(y + 1.0) + (x*x) + (y + 1.0) + 1.0
		
	}
	
	public static void test_setArgIdx() {
		MathFunc f = t + s;
		System.out.println(f);
		Map<String, MathFunc> fInners = new HashMap<String, MathFunc>();
		fInners.put("t", x*x);
		fInners.put("s", y*y);
		MathFunc fc = f.compose(fInners);
		check(fc.toString(), "f(x,y) = (x*x) + (y*y)");

		//Test different order of variable name list
		List<String> varNames = new ArrayList<String>();
		varNames.add("t");
		varNames.add("s");
		fc.setActiveVarNames(varNames);
		fc = fc + r;
		
		check(fc.toString(), "f(r,s,t) = t + s + r");
		System.out.println(fc.apply(new double[]{1.0,2.0,3.0})); //6.0
	}

	public static void test1() {
		MathFunc f = y - x;
		Map<String,Integer> argMap = new HashMap<String, Integer>();
		argMap.put("x", 1);
		argMap.put("y", 0);
		f.setArgIdx(argMap); //f(y,x)=y-x
		
		System.out.println(f.toString());
		System.out.println(f.getVarNames());
		System.out.println(f.getArgIdxMap());
		System.out.println(f.apply(1,3));

		MathFunc fx = r;
		MathFunc fy = s;
		Map<String, MathFunc> map = new HashMap<String, MathFunc>();
		map.put("x", fx); //x(r)=r
		map.put("y", fy); //y(s)=s
		MathFunc cf = f.compose(map);
		System.out.println(cf); //f(r,s)=s-r
		System.out.println(cf.getVarNames());
		System.out.println(cf.getArgIdxMap());
		
		//the order of "s","r" matters
		CompiledFunc ccf = cf.compile(new String[]{"s","r"});
		System.out.println(ccf.apply(5.5, 2));

	}
	
	public static void main(String[] args) {
		//test0();
		//test_setArgIdx();
		test1();
	}


}
