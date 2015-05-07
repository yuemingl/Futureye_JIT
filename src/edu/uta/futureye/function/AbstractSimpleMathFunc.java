package edu.uta.futureye.function;

import java.util.ArrayList;
import java.util.List;

import edu.uta.futureye.function.intf.MathFunc;

/**
 * Template to implement single variable MathFunc 
 *
 */
public abstract class AbstractSimpleMathFunc extends MathFuncBasic {
	protected String varName = null;
	protected int argIdx = 0;
	
	public AbstractSimpleMathFunc(String varName) {
		this.varName = varName;
	}

	@Override
	public String getName() {
		return varName;
	}

	@Override
	public MathFunc setName(String name) {
		this.varName = name;
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

	@Override
	public MathFunc setArgIdx(int... argIdx) {
		this.argIdx = argIdx[0];
		return this;
	}
	@Override
	public boolean isConstant() {
		return false;
	}
}
