package edu.uta.futureye.test;

import static edu.uta.futureye.function.FMath.*;
import static java.lang.Math.*;

public class TestMathematicalFunctions {
	public static double eps = 1e-5;
	public static void check(String info, double d1, double d2) {
		if(Math.abs(d1-d2) < eps) {
			System.out.println("pass");
		} else {
			System.out.println("!!!FAIL!!!   "+d1+"!="+d2+" "+info);
		}
	}
	
	public static void test1() {
		System.out.println(sin(x)+" "+sin(x).diff("x").setName("f'"));
		System.out.println(cos(x)+" "+cos(x).diff("x").setName("f'"));
		System.out.println(tan(x)+" "+tan(x).diff("x").setName("f'"));
		System.out.println(abs(x)+" "+abs(x).diff("x").setName("f'"));
		System.out.println(signum(x));
		System.out.println(sinh(x));
		System.out.println(cosh(x));
		System.out.println(tanh(x));
		System.out.println(asin(x));
		System.out.println(acos(x));
		System.out.println(log(x)+" "+log(x).diff("x").setName("f'"));
		System.out.println(log10(x)+" "+log10(x).diff("x").setName("f'"));
		System.out.println(max(x,y));
		System.out.println(min(x,y));
		System.out.println(exp(x)+" "+exp(x).diff("x").setName("f'"));
		System.out.println(pow(x,y)+" "+pow(x,y).diff("x").setName("f'"));
	}
	
	public static void test2() {
		check(sin(x).toString(), sin(x).apply(0.1), sin(0.1));
		check(cos(x).toString(), cos(x).apply(0.1), cos(0.1));
		check(tan(x).toString(), tan(x).apply(0.1), tan(0.1));
		check(abs(x).toString(), abs(x).apply(-0.1), abs(-0.1));
		check(signum(x).toString(), signum(x).apply(0.1), signum(0.1));
		check(sinh(x).toString(), sinh(x).apply(0.1), sinh(0.1));
		check(cosh(x).toString(), cosh(x).apply(0.1), cosh(0.1));
		check(asin(x).toString(), asin(x).apply(0.1), asin(0.1));
		check(acos(x).toString(), acos(x).apply(0.1), acos(0.1));
		check(log(x).toString(), log(x).apply(0.1), log(0.1));
		check(log10(x).toString(), log10(x).apply(0.1), log10(0.1));
		check(max(x,y).toString(), max(x,y).apply(0.1,0.2), max(0.1,0.2));
		check(min(x,y).toString(), min(x,y).apply(0.1,0.2), min(0.1,0.2));
		check(exp(x).toString(), exp(x).apply(0.1), exp(0.1));
		check(pow(x,y).toString(), pow(x,y).apply(0.1,0.2), pow(0.1,0.2));
	}
	
	public static void test3() {
		check(sin(x).toString(), sin(x).compile().apply(0.1), sin(0.1));
		check(cos(x).toString(), cos(x).compile().apply(0.1), cos(0.1));
		check(tan(x).toString(), tan(x).compile().apply(0.1), tan(0.1));
		check(abs(x).toString(), abs(x).compile().apply(0.1), abs(0.1));
		check(signum(x).toString(), signum(x).compile().apply(0.1), signum(0.1));
		check(sinh(x).toString(), sinh(x).compile().apply(0.1), sinh(0.1));
		check(cosh(x).toString(), cosh(x).compile().apply(0.1), cosh(0.1));
		check(asin(x).toString(), asin(x).compile().apply(0.1), asin(0.1));
		check(acos(x).toString(), acos(x).compile().apply(0.1), acos(0.1));
		check(log(x).toString(), log(x).compile().apply(0.1), log(0.1));
		check(log10(x).toString(), log10(x).compile().apply(0.1), log10(0.1));
		check(max(x,y).toString(), max(x,y).compile().apply(0.1,0.2), max(0.1,0.2));
		check(min(x,y).toString(), min(x,y).compile().apply(0.1,0.2), min(0.1,0.2));
		check(exp(x).toString(), exp(x).compile().apply(0.1), exp(0.1));
		check(pow(x,y).toString(), pow(x,y).compile().apply(0.1,0.2), pow(0.1,0.2));
	}	
	
	
	public static void main(String[] args) {
		test1();
		test2();
		test3();
	}
}
