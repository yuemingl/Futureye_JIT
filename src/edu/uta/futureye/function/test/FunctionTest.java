package edu.uta.futureye.function.test;

import edu.uta.futureye.util.Constant;

public class FunctionTest {
	
	public static void symbolTest() {
		System.out.println("\n begin symbol test...");

		FAxpb f1 = new FAxpb("f1(x)",2.0,4.0);
		FAxpb f2 = new FAxpb("f2(x)",1.0,2.0);
		System.out.println("f1="+f1);
		System.out.println("f2="+f2);
		
		System.out.println("2.0*f1+3.0*f2="+FOperator.LinearCombination(2.0, f1, 3.0, f2));
		System.out.println("f1^2="+FOperator.Power(f1, 2));
		System.out.println("3*f1="+FOperator.ScalarProduct(3.0, f1));

		System.out.println("f1+f2="+FOperator.Plus(f1, f2));
		System.out.println("f1-f2="+FOperator.Minus(f1, f2));
		System.out.println("f1*f2="+FOperator.Multi(f1, f2));
		System.out.println("f1/f2="+FOperator.Divi(f1, f2));
		
		System.out.println("end!");		
	}
	
	public static void valueTest() {
		System.out.println("\n begin value test...");

		FAxpb f1 = new FAxpb("F(x)=a*x+b",2.0,4.0);
		FAxpb f2 = new FAxpb("x",1.0,2.0);
		System.out.println("f1="+f1);
		System.out.println("f1(3.0)="+f1.getValue(new Variable("x",3.0)));

		System.out.println("f2="+f2);
		System.out.println("f2(3.0)="+f2.getValue(new Variable("x",3.0)));
		
		Function f1_dx = (Function) f1._d("x");
		System.out.println("f1_dx="+f1_dx);
		System.out.println("f1_dx(3.0)="+f1_dx.getValue(new Variable("x",3.0)));
		
		Function lc = FOperator.LinearCombination(2.0, f1, 3.0, f2);
		System.out.println("2.0*f1+3.0*f2="+lc);
		System.out.println("f(x=3.0)="+lc.getValue(new Variable("x",3.0)));
		
		Function power = FOperator.Power(f1, 2);
		System.out.println("f1^2="+power);
		System.out.println("f(x=3.0)="+power.getValue(new Variable("x",3.0)));
		
		Function sp = FOperator.ScalarProduct(3.0, f1);
		System.out.println("3*f1="+sp);
		System.out.println("f(x=3.0)="+sp.getValue(new Variable("x",3.0)));

		Function plus = FOperator.Plus(f1, f2);
		System.out.println("f1+f2="+plus);
		System.out.println("f(x=3.0)="+plus.getValue(new Variable("x",3.0)));
		
		Function minus = FOperator.Minus(f1, f2);
		System.out.println("f1-f2="+minus);
		System.out.println("f(x=3.0)="+minus.getValue(new Variable("x",3.0)));
		
		Function multi = FOperator.Multi(f1, f2);
		System.out.println("f1*f2="+multi);
		System.out.println("f(x=3.0)="+multi.getValue(new Variable("x",3.0)));
		
		Function divi = FOperator.Divi(f1, f2);
		System.out.println("f1/f2="+divi);
		System.out.println("f(x=3.0)="+divi.getValue(new Variable("x",3.0)));
		
		System.out.println("end!");		
	}
	
	public static void deritiveTest() {
		System.out.println("\n begin deritive test...");
		
		FAxpb f1 = new FAxpb("F1(x)=",5.0,7.0);
		FAxpb f2 = new FAxpb("F2(x)=",6.0,2.0);
		System.out.println("f1="+f1);
		System.out.println("f2="+f2);
		
		
		System.out.println("f_dx="+f1._d("x"));
		System.out.println("f_dx="+f2._d("x"));

		Function multi = FOperator.Multi(f1, f2);
		System.out.println("f1*f2="+multi);
		Function f_dx = (Function) multi._d("x");
		System.out.println("(f1*f2)_dx="+f_dx);
		System.out.println("(f1*f2)_dx="+f_dx.expand());
		
		System.out.println("end!");		
	}
	
	public static void chainOperatorTest() {
		System.out.println("\n begin chainOperator test...");
		
		FAxpb f1 = new FAxpb("F1(x)=",3.0,7.0);
		FAxpb f2 = new FAxpb("F2(x)=",5.0,11.0);
		Function fp = FOperator.Power(f2, 2.0);
		
		Function multi = FOperator.Multi(f1, f2);
		
		//multi = FOperator.Multi(multi, f2);
		//multi = FOperator.Multi(multi, f2);

		multi = FOperator.Multi(multi, fp);
		
		System.out.println("f1*f2="+multi);
		Function f12 = multi.expand();
		System.out.println("f1*f2="+f12);
		System.out.println("(f1*f2)_dx="+f12._d("x"));
		Function f12_d_e = (Function)((Function)f12._d("x")).expand();
		System.out.println("(f1*f2)_dx="+f12_d_e);
		
		System.out.println("(f1*f2)(1.0) = "+f12.getValue(new Variable("x",1.0)));
		
		Function divi = FOperator.Divi(f1, f2);
		System.out.println("f1*f2="+divi);
		Function f1d2 = divi.expand();
		System.out.println("f1/f2="+f1d2);
		System.out.println("(f1/f2)_dx="+f1d2._d("x"));
		System.out.println("(f1/f2)_dx="+((Function)f1d2._d("x")).expand());
		
		System.out.println("end!");
	}
	
	public static void composeTest() {
		System.out.println("\n begin compose test...");

		FAxpb f1 = new FAxpb("f1(x)",2.0,4.0);
		FAxpb f2 = new FAxpb("f2(r)","r",3.0,2.0);
		System.out.println("f1="+f1);
		System.out.println("f2="+f2);
		
		Function f12 = f1.compose(new ComposePair(f2,"x"));
		System.out.println("f12=f1(f2(x))="+f12);
		System.out.println("f1="+f1);
		System.out.println("f2="+f2);
		System.out.println("f12(3.0)="+f12.getValue(new Variable("r",3.0)));
		
		Function f12_dx = (Function)f12._d("x");
		System.out.println("f12_dx="+f12_dx);
		
		Function f12_dr = (Function)f12._d("r");
		System.out.println("f12_dr="+f12_dr);
		System.out.println("f12_dr="+f12_dr.expand());
		
		System.out.println("end!");
	}
	
	public static void TriFunction() {
		System.out.println("\n begin compose test...");
		
		FSin sin = new FSin("x");
		FCos cos = new FCos("x");
		System.out.println("f1(x) = "+sin);
		System.out.println("f1(x)_dx="+sin._d("x"));
		System.out.println("f2(x) = "+cos);
		System.out.println("f1(x)_dx="+cos._d("x"));

		Function sc = new AbstractFunction() {
			@Override
			public void createChain() {
				chain = new PlusChain();
				chain.addItem(new ItemPair(1.0,FOperator.Power(new FSin("x"),2.0)));
				chain.addItem(new ItemPair(1.0,FOperator.Power(new FCos("x"),2.0)));
			}
		};

		System.out.println("sc = "+sc);
		for(int i=1;i<10;i++) {
			System.out.println("sc("+Math.PI/i+")="+sc.getValue(new Variable("x",Math.PI/i)));
		}
		System.out.println("sc_dx = "+sc._d("x"));
		System.out.println("sc_dx = "+((Function)sc._d("x")).expand());
		
		FAxpb f = new FAxpb("f(r)","r",2.0,4.0);
		System.out.println("f(r)="+f);
		Function sinf = sin.compose(new ComposePair(f,"x"));
		System.out.println("sin(x(r))=" + sinf);
		Function sinfdr = (Function)sinf._d("r");
		System.out.println("sin(x(r))_dr=" + sinfdr.expand());
		
		System.out.println("end!");
		
	}
	
	public static void fullTest() {
		FX f = new FX("f1");
		for(int i=-1000;i<=1000;i++) {
			double v = Math.random()*i;
			double fv = f.getValue(new Variable("x",v));
			if(Math.abs(v-fv) > Constant.eps)
				System.out.println("Error: f(x)=x, x="+v+", f(x)="+fv);
		}
		
		
		
		
		
		
	}
	
	public static void serendipityShapFunction() {
		FX fx = new FX("r","r");
		FX fy = new FX("s","s");
		FAxpb f1mx = new FAxpb("1-r","r",-1.0,1.0);
		FAxpb f1px = new FAxpb("1+r","r",1.0,1.0);
		FAxpb f1my = new FAxpb("1-s","s",-1.0,1.0);
		FAxpb f1py = new FAxpb("1+s","s",1.0,1.0);

		Function N1 = FOperator.Multi(
			FOperator.Multi(FOperator.ScalarProduct(-0.25, f1mx), 
			f1my),
			FOperator.Plus(f1px, fy));
		
		
		System.out.println(N1);
		System.out.println(N1.expand().expand());
		//System.out.println(N1._d("r"));
		Function a = N1.expand().expand();
		System.out.println(((Function)a._d("r")).expand());
		System.out.println(a._d("s"));
	}	
	
	public static void main(String[] args) {
		symbolTest();
		valueTest();
		deritiveTest();
		chainOperatorTest();
		composeTest();
		TriFunction();
//		fullTest();
		serendipityShapFunction();
	}
}
