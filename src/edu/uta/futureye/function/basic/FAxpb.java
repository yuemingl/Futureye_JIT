package edu.uta.futureye.function.basic;

import java.util.Map;

import edu.uta.futureye.function.SingleVarFunc;
import edu.uta.futureye.function.FMath;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.VariableArray;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.util.Constant;

/**
 * f(x) = a*x + b
 */
public class FAxpb extends SingleVarFunc {
	protected double a;
	protected double b;

	public FAxpb(double a, double b) {
		super("", Constant.x);
		this.a = a;
		this.b = b;
	}
	
	public FAxpb(String varName, double a, double b) {
		super("", varName);
		this.a = a;
		this.b = b;
	}
	
	@Override
	public MathFunc diff(String varName) {
		if(this.varName.equals(varName))
			return new FC(a);
		else
			return FMath.C0;
	}

	@Override
	public double apply(Variable v) {
		return a*v.get(varName)+b;
	}

	@Override
	public double apply(double... args) {
		return a*args[argIdx]+b;
	}
	
	@Override
	public double apply(Variable v, Map<Object,Object> cache) {
		return a*v.get(varName)+b;
	}
	
	@Override
	public double[] applyAll(VariableArray v, Map<Object,Object> cache) {
		int len = v.length();
		double[] rlt = new double[len];
		double[] vs = v.get(varName);
		for(int i=0;i<len;i++)
			rlt[i] = a*vs[i]+b;
		return rlt;
	}
	
	@Override
	public int getOpOrder() {
		if(Double.compare(a, 0.0) == 0)
			return OP_ORDER0;
		if(Double.compare(b, 0.0) == 0)
			return OP_ORDER2;
		else
			return OP_ORDER3;
	}
	
	public String getExpr() {
		if(Double.compare(a, 1.0) == 0) {
			if(Double.compare(b, 0.0) == 0)
				return varName;
			else
				return varName+"+"+b;
		} else if(Double.compare(a, 0.0) == 0) {
				return b+"";
		} else if(Double.compare(b, 0.0) == 0) {
			return a+"*"+varName;
		}
		return a+"*"+varName+"+"+b;
	}

//TODO	
//	/**
//	 * varNames在由单个自变量表达式运算组合而成的多自变量函数情况下计算导数后可能被修改，
//	 * 这种修改是允许的，例如：(x+1)*(y+1)关于x求导数后，原来(y+1)的varNames由[y]变为[x,y]
//	 * 
//	 */
//	@Override
//	public MathFunc setVarNames(List<String> varNames) {
//		this.varNames = varNames;
//		return this;
//	}
}
