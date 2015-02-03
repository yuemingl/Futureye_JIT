package edu.uta.futureye.tutorial;

import static edu.uta.futureye.function.operator.FMath.*;

import java.util.HashMap;
import java.util.Map;

import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.FAxpb;
import edu.uta.futureye.function.intf.MathFun;

/**
 * Demo for basic function operations
 * 
 * @author liuyueming
 *
 */
public class T01Start {

	public static void main(String[] args) {
		//Definition of variables
		System.out.println("Definition of variables of functions: ");
		Variable v = new Variable(2.0);//x=2.0
		System.out.println(v);
		Variable v2 = new Variable("x",2.0).set("y", 3.0);//x=2.0 y=3.0
		System.out.println(v2);
		
		//Constant function
		System.out.println("\nConstant function: ");
		System.out.println("C(2.0)="+C(2.0));
		
		//f(x)=2.0*x+3.0
		System.out.println("\nDifferent ways of creating functions: ");
		MathFun f = C(2.0).M(X).A(C(3.0));
		System.out.println("f(x)="+f);
		//Evaluate function f(x) at variable v (x=2.0)
		System.out.println("Evaluate: f("+v.get()+")="+f.apply(v));
		
		//Another way to create f(x)=2.0*x+3.0
		MathFun f2 = new FAxpb(2.0,3.0);
		System.out.println("f(x)="+f2);
		System.out.println("Evaluate: f("+v.get()+")="+f2.apply(v));
		
		//Derivative of f(x): f'(x)=df(x)/dx
		System.out.println("\nDerivative of f(x):");
		System.out.println("df(x)/dx="+f._d("x"));
		
		//Composite functions
		System.out.println("\nComposite functions:");
		MathFun fOut = X.M(X).S(C1); //fOut(x)=x*x-1
		System.out.println("fOut(x)="+fOut);
		MathFun fIn = R.M(R).A(C1); //fIn(r)=r*r+1
		System.out.println("fIn(r)="+fIn);
		
		Map<String,MathFun> map = new HashMap<String,MathFun>();
		map.put("x", fIn);
		//fComp = x(r)*x(r) - 1
		MathFun fComp = fOut.compose(map);
		System.out.println("fComp( x(r) )="+fComp);
		System.out.println("Evaluate at x=2.0: fComp( 2.0 )="+
				fComp.apply(new Variable("x",2.0)));
		System.out.println("Evaluate at r=2.0: fComp( x(2.0) )="+
				fComp.apply(new Variable("r",2.0)));
		
		System.out.println("d(fComp)/dx="+fComp._d("x"));
		MathFun dr = fComp._d("r");
		System.out.println("d(fComp)/dr="+dr);
		System.out.println("d(fComp)/dr|_r=2.0  = "+
				dr.apply(new Variable("r",2.0)));
	}

}
