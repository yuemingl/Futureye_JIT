package edu.uta.futureye.function.basic;

import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.function.AbstractSimpleMathFunc;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.util.Constant;

/**
 * f(x) = a*x
 *
 */
public class FAx extends AbstractSimpleMathFunc {
	protected double a;
	
	public FAx(double a) {
		super(Constant.x);
		varName = Constant.x;
		this.a = a;
	}
	
	public FAx(String varName, double a) {
		super(varName);
		this.a = a;
	}
	
	@Override
	public MathFunc diff(String varName) {
		if(this.getVarNames().contains(varName))
			return new FC(a);
		else
			return FC.C0;
	}

	@Override
	public double apply(Variable v) {
		return a*v.get(getVarNames().get(0));
	}

	@Override
	public double apply(Element e, Node n, double... args) {
		return a*args[0];
	}

	@Override
	public double apply(double... args) {
		return apply(null, null, args);
	}
	
	@Override
	public MathFunc copy() {
		return new FAx(this.varName, a);
	}

	@Override
	public String getExpr() {
		return toString();
	}
	
	@Override
	public String toString() {
		if(Double.compare(a, 1.0) == 0)
			return " "+getVarNames().get(0)+" ";
		else if(Double.compare(a, 0.0) == 0)
			return " 0.0 ";
		return " "+a+"*"+getVarNames().get(0)+" ";
	}

	@Override
	public boolean isConstant() {
		return false;
	}
}
