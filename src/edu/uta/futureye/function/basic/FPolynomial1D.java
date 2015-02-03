package edu.uta.futureye.function.basic;

import java.util.List;

import edu.uta.futureye.function.AbstractMathFun;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.intf.MathFun;
import edu.uta.futureye.util.Constant;

/**
 * f(x) = an*x^n + an_1*x^(n-1) + ... + a1*x + a0
 * @author liuyueming
 *
 */
public class FPolynomial1D extends AbstractMathFun {
	List<Double> coefList;
	
	/**
	 * 构造一个多项式
	 * @param coefList
	 * a0 = coefList.get(0)
	 * ...
	 * an = coefList.get(coefList.size()-1)
	 */
	public FPolynomial1D(List<Double> coefList) {
		varNames.add(Constant.x);
		this.coefList = coefList;
	}
	
	@Override
	public MathFun _d(String varName) {
		if(this.getVarNames().contains(varName))
			return derivative1(1,1);
		else 
			return new FC(0.0);
	}
	
	protected FPolynomial1D derivative1(int degree,int maxDegree) {
		int d = degree;
		if(d >= coefList.size()) {
			coefList.clear();
			return this;
		}
		if(0<d && d <= maxDegree) {
			coefList.remove(0);
			for(int i=0;i<coefList.size();i++) {
				coefList.set(i, coefList.get(i)*(i+1));
			}
			derivative1(--d,maxDegree);
		}
		return this;
	}
	
	@Override
	public double apply(Variable v) {
		double x = v.get(getVarNames().get(0));
		double f = 0.0;
		for(int i=0;i<coefList.size();i++) {
			f += coefList.get(i)*Math.pow(x, i);
		}
		return f;
	}
}
