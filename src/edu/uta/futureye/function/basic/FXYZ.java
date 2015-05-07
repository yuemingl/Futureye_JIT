package edu.uta.futureye.function.basic;

import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.function.AbstractMathFun;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.util.Constant;

/**
 * f(x,y) = c1*x + c2*y + c3*z + c4
 * 
 */
public class FXYZ extends AbstractMathFun{

	public FXYZ() {
		varNames.add(Constant.x);
		varNames.add(Constant.y);
		varNames.add(Constant.z);
	}
	
	public FXYZ(String varName) {
		super(varName);
	}

	
	@Override
	public MathFunc diff(String varName) {
		return null;
	}

	@Override
	public double apply(Variable v) {
		return v.get(getVarNames().get(0));
	}
	
	public String toString() {
		return getVarNames().get(0);
	}

	@Override
	public double apply(Element e, Node n, double... args) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public double apply(double... args) {
		return apply(null, null, args);
	}
}