package edu.uta.futureye.function;

import java.util.LinkedList;
import java.util.List;

import edu.uta.futureye.function.intf.MathFunc;

/**
 * Template to implement multiple variables MathFunc 
 *
 */
public abstract class AbstractMathFunc extends MathFuncBasic {
	protected List<String> varNames = new LinkedList<String>();
	protected int[] argIdx;
	protected String fName = null;
	
	public AbstractMathFunc() {
	}
	
	public AbstractMathFunc(List<String> varNames) {
		this.varNames = varNames;
	}
	
	public AbstractMathFunc(String varName, String ...aryVarNames) {
		varNames.add(varName);
		for(String s : aryVarNames)
			varNames.add(s);
	}
	
	@Override
	public List<String> getVarNames() {
		return varNames;
	}

	@Override
	public MathFunc setVarNames(List<String> varNames) {
		this.varNames = varNames;
		return this;
	}

	public MathFunc setArgIdx(int ...argIdx) {
		this.argIdx = argIdx;
		return this;
	}

	@Override
	public String getName() {
		return this.fName;
	}
	
	@Override
	public MathFunc setName(String name) {
		this.fName = name;
		return this;
	}
	
	@Override
	public boolean isConstant() {
		return false;
	}
	
	@Override
	public String getExpr() {
		String varList = getVarNames().toString();
		String displayVarList = "("+varList.substring(1, varList.length()-1)+")";
		
		Class<?> enclosingClass = getClass().getEnclosingClass();
		if (enclosingClass != null) {
		  return enclosingClass.getSimpleName() + displayVarList;
		} else {
		  return getClass().getSimpleName() + displayVarList;
		}
	}
	
	@Override
	public String toString() {
		if(getName() == null) {
			return getExpr();
		} else 
			return getName();
	}
}
