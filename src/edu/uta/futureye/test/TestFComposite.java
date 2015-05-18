package edu.uta.futureye.test;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import edu.uta.futureye.function.basic.FX;
import edu.uta.futureye.function.intf.MathFunc;
import static edu.uta.futureye.function.FMath.*;

public class TestFComposite {
	public static void test1() {
		MathFunc f = FX.r.A(FX.s);
		System.out.println(f);
		Map<String, MathFunc> fInners = new HashMap<String, MathFunc>();
		fInners.put("r", FX.x*FX.x);
		fInners.put("s", FX.y*FX.y);
		MathFunc fc = f.compose(fInners);
		System.out.println(fc);
		
		double[] args = new double[]{2.0, 3.0};
		System.out.println(fc.apply(args));
		System.out.println(fc.compile().apply(args));
		
		List<String> varNames = new ArrayList<String>();
		varNames.add("r");
		varNames.add("s");
		
		fc.setActiveVarNames(varNames);
		System.out.println(fc.apply(args));
		System.out.println(fc.compile().apply(args));
	}
	
	public static void test2() {
		MathFunc f = r*s + 1;
		Map<String, MathFunc> fInners = new HashMap<String, MathFunc>();
		fInners.put("r", x*x);
		fInners.put("s", y+1);
		MathFunc fc = f.compose(fInners);
		System.out.println(fc); //
	}

	public static void main(String[] args) {
		test1();
		test2();
	}


}
