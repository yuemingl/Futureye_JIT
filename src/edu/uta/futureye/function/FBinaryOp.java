package edu.uta.futureye.function;

import java.util.Map;

import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.util.Utils;

public abstract class FBinaryOp extends AbstractMathFunc {
	public MathFunc arg1;
	public MathFunc arg2;
	
	public FBinaryOp(MathFunc arg1, MathFunc arg2) {
		this.arg1 = arg1;
		this.arg2 = arg2;
	}
	
	public MathFunc[] args() {
		return new MathFunc[] { arg1, arg2 };
	}
	
	public MathFunc lhs() {
		return arg1;
	}
	
	public MathFunc rhs() {
		return arg2;
	}
	
	@Override
	public MathFunc setArgIdx(Map<String, Integer> argsMap) {
		for(int i=0; i<varNames.length; i++) {
			this.argIdx[i] = argsMap.get(varNames[i]);
		}
		
		if(!Utils.isContained(argsMap, this.arg1.getArgIdxMap()))
			this.arg1 = this.arg1.copy().setArgIdx(argsMap);
		
		if(!Utils.isContained(argsMap, this.arg2.getArgIdxMap()))
			this.arg2 = this.arg2.copy().setArgIdx(argsMap);
		
		return this;
	}
}
