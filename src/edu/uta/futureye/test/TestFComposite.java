package edu.uta.futureye.test;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import edu.uta.futureye.function.basic.FX;
import edu.uta.futureye.function.intf.MathFunc;
import static edu.uta.futureye.function.FMath.*;

public class TestFComposite {
	public static void test0() {
		MathFunc f = r*s + r + s + 1;
		Map<String, MathFunc> fInners = new HashMap<String, MathFunc>();
		fInners.put("r", x*x);
		fInners.put("s", y+1);
		MathFunc fc = f.compose(fInners);
		System.out.println(fc); //f(x,y) = (x*x)*(y + 1.0) + (x*x) + (y + 1.0) + 1.0
		System.out.println(fc.compile().apply(new double[]{2.0,3.0})); //25.0
		System.out.println(fc.compile(new String[]{"a","b","x","y"}).apply(new double[]{1.0,2.0,2.0,3.0}));
	}
	
	public static void test1() {
		MathFunc f = r + s;
		System.out.println(f);
		Map<String, MathFunc> fInners = new HashMap<String, MathFunc>();
		fInners.put("r", x*x);
		fInners.put("s", y*y);
		MathFunc fc = f.compose(fInners);
		System.out.println(fc);
		
		double[] args = new double[]{2.0, 3.0};
		System.out.println(fc.apply(args));
		System.out.println(fc.compile().apply(args)); //13.0
		System.out.println(fc.compile(new String[]{"a","b","x","y"}).apply(new double[]{1.0,2.0,2.0,3.0}));

	}
	
	public static void test_setArgIdx() {
		MathFunc f = t + s;
		System.out.println(f);
		Map<String, MathFunc> fInners = new HashMap<String, MathFunc>();
		fInners.put("t", FX.x*FX.x);
		fInners.put("s", FX.y*FX.y);
		MathFunc fc = f.compose(fInners);
		System.out.println(fc);

		//Test different order of variable name list
		List<String> varNames = new ArrayList<String>();
		varNames.add("t");
		varNames.add("s");
		fc.setActiveVarNames(varNames);
		fc = fc + r;
		
		System.out.println(fc);
		System.out.println(fc.apply(new double[]{1.0,2.0,3.0})); //6.0
	}


	public static void main(String[] args) {
		test0();
		test1();
		test_setArgIdx();
	}


}
