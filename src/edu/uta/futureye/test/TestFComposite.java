package edu.uta.futureye.test;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

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

	public static void main(String[] args) {
		test0();
		test_setArgIdx();
	}


}
