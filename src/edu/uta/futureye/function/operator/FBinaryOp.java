package edu.uta.futureye.function.operator;

import java.util.List;
import java.util.Map;

import edu.uta.futureye.function.MultiVarFunc;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.util.Utils;

public abstract class FBinaryOp extends MultiVarFunc {
	public MathFunc arg1;
	public MathFunc arg2;
	
	public FBinaryOp(MathFunc arg1, MathFunc arg2) {
		this.arg1 = arg1;
		this.arg2 = arg2;
		List<String> list = Utils.mergeList(arg1.getVarNames(), arg2.getVarNames());
		Map<String, Integer> map = Utils.getIndexMap(list);
		setVarNames(list);
		setArgIdx(map);
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
	
	/**
	 * TODO return this if no change
	 * return the copied version if there is any change for argIdx
	 * this.arg1 = this.arg1.setArgIdx(...)
	 * this.arg2 = this.arg2.setArgIdx(...)
	 * 
	 */
	@Override
	public MathFunc setArgIdx(Map<String, Integer> argsMap) {
		//Allocate new array each time due to the "copy on change"
		int[] idx = new int[varNames.length];
		for(int i=0; i<varNames.length; i++) {
			idx[i] = argsMap.get(varNames[i]);
		}
		this.argIdx = idx;
		
		//Copy on change
		if(!Utils.isMapContain(argsMap, this.arg1.getArgIdxMap()))
			this.arg1 = this.arg1.copy().setArgIdx(argsMap);
		
		if(!Utils.isMapContain(argsMap, this.arg2.getArgIdxMap()))
			this.arg2 = this.arg2.copy().setArgIdx(argsMap);
		
		return this;
	}
}
