package edu.uta.futureye.function.basic;

import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.function.SingleVarFunc;
import edu.uta.futureye.function.FMath;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.util.Constant;

/**
 * f(x) = a*x
 *
 */
public class FAx extends SingleVarFunc {
	protected double a;
	
	public FAx(double a) {
		super(null, Constant.x);
		this.a = a;
	}
	
	public FAx(String varName, double a) {
		super("", varName);
		this.a = a;
	}
	
	@Override
	public MathFunc diff(String varName) {
		if(this.getVarNames().contains(varName))
			return new FC(a);
		else
			return FMath.C0;
	}

	@Override
	public double apply(Variable v) {
		return a*v.get(getVarNames().get(0));
	}

	@Override
	public double apply(Element e, Node n, double... args) {
		return a*args[argIdx];
	}

	@Override
	public double apply(double... args) {
		return a*args[argIdx];
	}
	
//	@Override
//	public MathFunc copy() {
//		FAx ret = new FAx(this.varName, a);
//		ret.fName = this.fName;
//		ret.argIdx = this.argIdx;
//		return ret;
//	}

	@Override
	public String getExpr() {
		if(Double.compare(a, 1.0) == 0)
			return varName;
		else if(Double.compare(a, 0.0) == 0)
			return "0.0";
		return a+"*"+varName;
	}
}
