package edu.uta.futureye.function.basic;

import java.util.List;

import edu.uta.futureye.function.MultiVarFunc;
import edu.uta.futureye.function.FMath;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.util.Constant;

/**
 * f(x,y,z) = c1*x + c2*y + c3*z + c4
 * 
 */
public class FXYZ extends MultiVarFunc{
	protected double c1, c2, c3, c4;

	public FXYZ(double c1, double c2, double c3, double c4) {
		this.varNames = new String[3];
		this.varNames[0] = Constant.x;
		this.varNames[1] = Constant.y;
		this.varNames[2] = Constant.z;
		this.c1 = c1;
		this.c2 = c2;
		this.c3 = c3;
		this.c4 = c4;
	}
	
	/**
	 * Construct the function with variable names other than the default x,y,z
	 * @param funcName
	 * @param varNames
	 * @param c1
	 * @param c2
	 * @param c3
	 * @param c4
	 */
	public FXYZ(String funcName, List<String> varNames, double c1, double c2, double c3, double c4) {
		super(funcName, varNames);
		this.c1 = c1;
		this.c2 = c2;
		this.c3 = c3;
		this.c4 = c4;
	}

	@Override
	public double apply(double... args) {
		return c1*args[0] + c2*args[1] + c3*args[3] + c4;
	}

	@Override
	public MathFunc diff(String varName) {
		if(varNames[0].equals(varName))
			return new FC(c1);
		else if(varNames[1].equals(varName)) {
			return new FC(c2);
		} else if(varNames[2].equals(varName)) {
			return new FC(c3);
		}
		return FMath.C0;
	}

	public String toString() {
		String s1 = "";
		if(Math.abs(c1) > Constant.eps)
			s1 = s1 + "*";
		String s2 = "";
		if(Math.abs(c2) > Constant.eps)
			s2 = s2 + "*";
		if(s1.isEmpty() && !s2.isEmpty())
			return s2+this.varNames[1];
		else if(!s1.isEmpty() && s2.isEmpty())
			return s1+this.varNames[0];
		else
			return s1+this.varNames[0]+" + " +s2+this.varNames[1];
		//TODO c3
	}
}