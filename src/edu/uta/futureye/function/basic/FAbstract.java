package edu.uta.futureye.function.basic;

import java.util.LinkedList;
import java.util.List;

import edu.uta.futureye.function.DerivativeIndicator;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.intf.FunctionDerivable;

public abstract class FAbstract implements FunctionDerivable {
	private List<String> varNames = null;

	public FAbstract() {
		varNames = new LinkedList<String>();
		varNames.add("x");
	}
	
	public FAbstract(String varName, String ...aryVarNames) {
		varNames = new LinkedList<String>();
		varNames.add(varName);
		if(aryVarNames != null) {
			for(String s : aryVarNames)
				varNames.add(s);
		}
	}
	
	public FAbstract(List<String> varNames) {
		this.varNames = varNames;
	}
	
	@Override
	public FunctionDerivable derivative(DerivativeIndicator di) {
		return null;
	}
	
	@Override
	public double value(Variable v) {
		return Double.NaN;
	}

	@Override
	public List<String> varNames() {
		return varNames;
	}

	@Override
	public void setVarNames(List<String> varNames) {
		this.varNames = varNames;
	}

	public String toString() {
		String s = varNames.toString();
		return "F("+s.substring(1, s.length()-1)+")";
	}
}
