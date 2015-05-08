package edu.uta.futureye.function;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import edu.uta.futureye.function.intf.MathFunc;

/**
 * Template to implement single variable MathFunc 
 *
 */
public abstract class AbstractSimpleMathFunc extends MathFuncBasic {
	protected String varName;
	protected int argIdx;
	protected String fName;
	
	public AbstractSimpleMathFunc(String funcName, String varName) {
		this.fName = funcName;
		this.varName = varName;
	}

	@Override
	public String getName() {
		return fName;
	}

	@Override
	public MathFunc setName(String name) {
		this.fName = name;
		return this;
	}

	@Override
	public List<String> getVarNames() {
		List<String> ret = new ArrayList<String>();
		ret.add(varName);
		return ret;
	}

	@Override
	public MathFunc setVarNames(List<String> varNames) {
		this.varName = varNames.get(0);
		return this;
	}
	
	public String getVarName() {
		return this.varName;
	}
	
	public MathFunc setVarName(String varName) {
		this.varName = varName;
		return this;
	}

	@Override
	public MathFunc setArgIdx(Map<String, Integer> argsMap) {
		this.argIdx = argsMap.get(varName);
		return this;
	}
	
	@Override
	public Map<String, Integer> getArgIdxMap() {
		Map<String, Integer> ret = new HashMap<String, Integer>();
		ret.put(varName, argIdx);
		return ret;
	}
	
	@Override
	public boolean isConstant() {
		return false;
	}
	
	@Override
	public String getExpr() {
		return fName + "(" + varName + ")";
	}
}
