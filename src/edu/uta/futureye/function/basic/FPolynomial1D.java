package edu.uta.futureye.function.basic;

import java.util.List;

import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.function.AbstractMathFunc;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.util.Constant;

/**
 * f(x) = an*x^n + an_1*x^(n-1) + ... + a1*x + a0
 *
 */
public class FPolynomial1D extends AbstractMathFunc {
	List<Double> coefList;
	
	/**
	 * 构造一个多项式
	 * @param coefList
	 * a0 = coefList.get(0)
	 * ...
	 * an = coefList.get(coefList.size()-1)
	 */
	public FPolynomial1D(List<Double> coefList) {
		super(Constant.x);
		this.coefList = coefList;
	}
	
	@Override
	public MathFunc diff(String varName) {
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

	@Override
	public double apply(Element e, Node n, double... args) {
		double x = args[0];
		double f = 0.0;
		for(int i=0;i<coefList.size();i++) {
			f += coefList.get(i)*Math.pow(x, i);
		}
		return f;
	}

	@Override
	public double apply(double... args) {
		return apply(null, null, args);
	}
}
