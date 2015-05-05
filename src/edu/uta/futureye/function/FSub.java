package edu.uta.futureye.function;

import java.util.Map;

import com.sun.org.apache.bcel.internal.generic.ConstantPoolGen;
import com.sun.org.apache.bcel.internal.generic.DSUB;
import com.sun.org.apache.bcel.internal.generic.InstructionConstants;
import com.sun.org.apache.bcel.internal.generic.InstructionFactory;
import com.sun.org.apache.bcel.internal.generic.InstructionHandle;
import com.sun.org.apache.bcel.internal.generic.InstructionList;
import com.sun.org.apache.bcel.internal.generic.MethodGen;

import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.util.Constant;
import edu.uta.futureye.util.Utils;

public class FSub extends FBinaryOp {
	public FSub(MathFunc left, MathFunc right) {
		super(left, right);
		setVarNames(Utils.mergeList(left.getVarNames(), right.getVarNames()));
	}

	@Override
	public double apply(Variable v) {
		return arg1.apply(v) - arg2.apply(v);
	}
	
	@Override
	public double apply(Variable v, Map<Object,Object> cache) {
		return arg1.apply(v,cache) - arg2.apply(v,cache);
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
	public MathFunc _d(String varName) {
		return arg1._d(varName).S(arg2._d(varName)).setVarNames(this.getVarNames());
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
	public InstructionHandle bytecodeGen(MethodGen mg, ConstantPoolGen cp,
			InstructionFactory factory, InstructionList il,
			Map<String, Integer> argsMap, int argsStartPos) {
		arg1.bytecodeGen(mg, cp, factory, il, argsMap, argsStartPos);
		arg2.bytecodeGen(mg, cp, factory, il, argsMap, argsStartPos);
		return il.append(InstructionConstants.DSUB);
	}

}
