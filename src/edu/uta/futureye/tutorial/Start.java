package edu.uta.futureye.tutorial;

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
		
		//f(x)=2*x+3
		Function f = FC.c(2.0).X(FX.fx).P(FC.c(3.0));
		System.out.println("f(x)="+f);
		System.out.println("f("+v.get()+")="+f.value(v));
		
		//Another way to get f(x)=2*x+3
		Function f2 = new FAxpb(2.0,3.0);
		System.out.println("f(x)="+f2);
		System.out.println("f("+v.get()+")="+f2.value(v));
		
		//Derivative: f'(x)=df(x)/dx
		System.out.println("df(x)/dx="+f.d("x"));
		
	}

}
