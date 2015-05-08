package edu.uta.futureye.function;

import java.util.List;
import java.util.Map;

import com.sun.org.apache.bcel.internal.generic.ConstantPoolGen;
import com.sun.org.apache.bcel.internal.generic.InstructionConstants;
import com.sun.org.apache.bcel.internal.generic.InstructionFactory;
import com.sun.org.apache.bcel.internal.generic.InstructionHandle;
import com.sun.org.apache.bcel.internal.generic.InstructionList;
import com.sun.org.apache.bcel.internal.generic.MethodGen;

import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.util.Constant;
import edu.uta.futureye.util.Utils;

public class FSub extends FBinaryOp {
	public FSub(MathFunc left, MathFunc right) {
		super(left, right);
		List<String> list = Utils.mergeList(left.getVarNames(), right.getVarNames());
		Map<String, Integer> map = Utils.getIndexMap(list);
		setVarNames(list);
		setArgIdx(map);
	}

	@Override
	public double apply(Variable v) {
		return arg1.apply(v) - arg2.apply(v);
	}

	@Override
	public double apply(double... args) {
		return arg1.apply(null, null, args) - arg2.apply(null, null, args);
	}
	
	@Override
	public double apply(Variable v, Map<Object,Object> cache) {
		return arg1.apply(v,cache) - arg2.apply(v,cache);
	}
	
	@Override
	public double apply(Element e, Node n, double... args) {
		return arg1.apply(e,n,args) - arg2.apply(e, n, args);
	}
	
	@Override
	public double[] applyAll(VariableArray v, Map<Object,Object> cache) {
		int len = v.length();
		double[] la = arg1.applyAll(v,cache);
		double[] ra = arg2.applyAll(v,cache);
		for(int i=0;i<len;i++) {
			la[i] -= ra[i];
		}
		return la;
	}
	
	@Override
	public MathFunc diff(String varName) {
		//return arg1.diff(varName).S(arg2.diff(varName)).setVarNames(this.getVarNames());
		return arg1.diff(varName).S(arg2.diff(varName));
	}
	
	@Override
	public int getOpOrder() {
		return OP_ORDER3;
	}
	
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		if(! (arg1.isConstant() && Math.abs(arg1.apply()) < Constant.eps) ) {
			sb.append(arg1.toString());
		}
		sb.append(" - ");
		if(arg2.getOpOrder() >= OP_ORDER3)
			sb.append("(").append(arg2.toString()).append(")");
		else
			sb.append(arg2.toString());
		return sb.toString();
	}

	@Override
	public InstructionHandle bytecodeGen(String clsName, MethodGen mg,
			ConstantPoolGen cp, InstructionFactory factory,
			InstructionList il, Map<String, Integer> argsMap, int argsStartPos, 
			Map<MathFunc, Integer> funcRefsMap) {
		arg1.bytecodeGen(clsName, mg, cp, factory, il, argsMap, argsStartPos, funcRefsMap);
		arg2.bytecodeGen(clsName, mg, cp, factory, il, argsMap, argsStartPos, funcRefsMap);
		return il.append(InstructionConstants.DSUB);
	}
	
//	@Override
//	public MathFunc copy() {
//		return new FSub(this.arg1, this.arg2).setName(this.getName());
//	}
}
