package edu.uta.futureye.function;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import edu.uta.futureye.function.intf.MathFunc;

/**
 * Template to implement multiple variables MathFunc 
 *
 */
public abstract class AbstractMathFunc extends MathFuncBasic {
	protected String[] varNames;
	protected int[] argIdx;
	protected String fName;
	
	public AbstractMathFunc() {
	}
	
	public AbstractMathFunc(List<String> varNames) {
		this.varNames = varNames.toArray(new String[0]);
		this.argIdx = new int[this.varNames.length];
	}
	
	public AbstractMathFunc(String varName, String ...aryVarNames) {
		List<String> list = new ArrayList<String>();
		list.add(varName);
		for(String s : aryVarNames)
			list.add(s);
		this.varNames = list.toArray(new String[0]);
		this.argIdx = new int[this.varNames.length];
	}
	
	@Override
	public List<String> getVarNames() {
		List<String> list = new ArrayList<String>();
		for(String s : varNames)
			list.add(s);
		return list;
	}

	@Override
	public MathFunc setVarNames(List<String> varNames) {
		this.varNames = varNames.toArray(new String[0]);
		this.argIdx = new int[this.varNames.length];
		return this;
	}
	
	@Override
	public MathFunc setArgIdx(Map<String, Integer> argsMap) {
		for(int i=0; i<varNames.length; i++) {
			this.argIdx[i] = argsMap.get(varNames[i]);
		}
		return this;
	}
	
	@Override
	public Map<String, Integer> getArgIdxMap() {
		Map<String, Integer> ret = new HashMap<String, Integer>();
		for(int i=0; i<varNames.length; i++) {
			ret.put(varNames[i], argIdx[i]);
		}
		return ret;
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
		
//		Class<?> enclosingClass = getClass().getEnclosingClass();
//		if (enclosingClass != null) {
//			return enclosingClass.getSimpleName() + displayVarList;
//		} else {
//			return getClass().getSimpleName() + displayVarList;
//		}
		return getName() + displayVarList;
	}
}
