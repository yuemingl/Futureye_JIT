package edu.uta.futureye.tutorial;

import java.util.HashMap;
import java.util.Map;

import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.FAxpb;
import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.basic.FX;
import edu.uta.futureye.function.intf.Function;

public class Start {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		//x=2.0
		Variable v = new Variable(2.0);
		System.out.println(v);
		
		//Constant function
		System.out.println(FC.c(2.0));
		
		//f(x)=2.0*x+3.0
		Function f = FC.c(2.0).M(FX.fx).A(FC.c(3.0));
		System.out.println("f(x)="+f);
		//evaluate function at v (x=2.0)
		System.out.println("f("+v.get()+")="+f.value(v));
		
		//Another way to get f(x)=2.0*x+3.0
		Function f2 = new FAxpb(2.0,3.0);
		System.out.println("f(x)="+f2);
		System.out.println("f("+v.get()+")="+f2.value(v));
		
		//Derivative: f'(x)=df(x)/dx
		System.out.println("df(x)/dx="+f._d("x"));
		
		//Composite functions
		Function fx = FX.fx; //Use predefined static object
		Function fr = new FX("r");
		
		Function fOut = fx.M(fx).S(FC.c1);
		System.out.println(fOut);
		Function fIn = fr.M(fr).A(FC.c1);
		System.out.println(fIn);
		
		Map<String,Function> map = new HashMap<String,Function>();
		map.put("x", fIn);
		Function fComp = fOut.compose(map);
		System.out.println(fComp);
		System.out.println(fComp.value(new Variable("x",2.0)));
		System.out.println(fComp.value(new Variable("r",2.0)));
		
		System.out.println(fComp._d("x"));
		Function dr = fComp._d("r");
		System.out.println(dr);
		System.out.println(dr.value(new Variable("r",2.0)));
	}

}
