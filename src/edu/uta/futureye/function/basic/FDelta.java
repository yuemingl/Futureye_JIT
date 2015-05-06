package edu.uta.futureye.function.basic;

import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.function.AbstractMathFun;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.intf.MathFunc;

public class FDelta extends AbstractMathFun {
	Variable x0;
	double eps;
	double amp;

	/**
	 * Delta(x) = amp * e^(-0.25*|x-x0|^2/eps ) / (2*sqrt(PI*eps))
	 * 
	 * e.g. Function det = new FDelta(new Variable("x",89).set("y",76).set("z",191),1e-4,1e5);
	 * 
	 * @param x0
	 * @param eps
	 * @param amp
	 */
	public FDelta(Variable x0,double eps,double amp) {
		if(x0.getNameValuePairs().size() == 1) {
			varNames.add("x");
		} else if(x0.getNameValuePairs().size() == 2) {
			varNames.add("x");
			varNames.add("y");
		} else if(x0.getNameValuePairs().size() == 3) {
			varNames.add("x");
			varNames.add("y");		
			varNames.add("z");		
		}
		this.x0 = x0;
		this.eps = eps;
		this.amp = amp;
	}
	
	@Override
	public double apply(Variable x) {
		double d2 = 0.0;
		if(x0.getNameValuePairs().size() == 1) {
			double dx = x.get("x")-x0.get("x");
			d2 = dx*dx;
		} else if(x0.getNameValuePairs().size() == 2) {
			double dx = x.get("x")-x0.get("x");
			double dy = x.get("y")-x0.get("y");
			d2 = dx*dx+dy*dy;
		} else if(x0.getNameValuePairs().size() == 3) {
			double dx = x.get("x")-x0.get("x");
			double dy = x.get("y")-x0.get("y");
			double dz = x.get("z")-x0.get("z");
			d2 = dx*dx+dy*dy+dz*dz;
		}
		return amp*Math.exp(-d2/eps/4.0) / (2*Math.sqrt(Math.PI*eps));
	}

	@Override
	public double apply(Element e, Node n, double... args) {
		double d2 = 0.0;
		if(x0.getNameValuePairs().size() == 1) {
			double dx = args[0]-x0.get("x");
			d2 = dx*dx;
		} else if(x0.getNameValuePairs().size() == 2) {
			double dx = args[0]-x0.get("x");
			double dy = args[1]-x0.get("y");
			d2 = dx*dx+dy*dy;
		} else if(x0.getNameValuePairs().size() == 3) {
			double dx = args[0]-x0.get("x");
			double dy = args[1]-x0.get("y");
			double dz = args[2]-x0.get("z");
			d2 = dx*dx+dy*dy+dz*dz;
		}
		return amp*Math.exp(-d2/eps/4.0) / (2*Math.sqrt(Math.PI*eps));
	}
	
	public double apply() {
		return 0.0;
	}
	
	@Override
	public MathFunc _d(String varName) {
		return null;
	}
	
	public String toString() {
		return "Delta";
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		Variable x0 = new Variable();
		x0.set("x", 0.5);
		x0.set("y", 0.5);
		
		Variable x = new Variable();
		x.set("x", 0.6);
		x.set("y", 0.6);
		
		FDelta delta = new FDelta(x0,1e-4,1e5);
		System.out.println(delta.apply(x));
	}
}
