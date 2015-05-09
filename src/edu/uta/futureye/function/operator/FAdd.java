package edu.uta.futureye.function.operator;

import java.util.Map;

import com.sun.org.apache.bcel.internal.generic.ConstantPoolGen;
import com.sun.org.apache.bcel.internal.generic.InstructionConstants;
import com.sun.org.apache.bcel.internal.generic.InstructionFactory;
import com.sun.org.apache.bcel.internal.generic.InstructionHandle;
import com.sun.org.apache.bcel.internal.generic.InstructionList;
import com.sun.org.apache.bcel.internal.generic.MethodGen;

import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.VariableArray;
import edu.uta.futureye.function.intf.MathFunc;

public class FAdd extends FBinaryOp {
	public FAdd(MathFunc left, MathFunc right) {
		super(left, right);
	}

	@Override
	public double apply(double... args) {
		return arg1.apply(null, null, args) + arg2.apply(null, null, args);
	}

	@Override
	public double apply(Element e, Node n, double... args) {
		return arg1.apply(e,n,args) + arg2.apply(e, n, args);
	}
	
	@Override
	public double[] applyAll(VariableArray v, Map<Object,Object> cache) {
		int len = v.length();
		double[] la = arg1.applyAll(v,cache);
		double[] ra = arg2.applyAll(v,cache);
		for(int i=0;i<len;i++) {
			la[i] += ra[i];
		}
		return la;
	}

	@Override
	public MathFunc diff(String varName) {
		//return arg1.diff(varName).A(arg2.diff(varName)).setVarNames(this.getVarNames());
		return arg1.diff(varName).A(arg2.diff(varName));
	}
	
	@Override
	public int getOpOrder() {
		return OP_ORDER3;
	}
	
	@Override
	public String getExpr() {
		StringBuilder sb = new StringBuilder();
		sb.append(arg1.getExpr());
		sb.append(" + ");
		sb.append(arg2.getExpr());
		return sb.toString();
	}

	@Override
	public InstructionHandle bytecodeGen(String clsName, MethodGen mg,
			ConstantPoolGen cp, InstructionFactory factory,
			InstructionList il, Map<String, Integer> argsMap, int argsStartPos, 
			Map<MathFunc, Integer> funcRefsMap) {
		arg1.bytecodeGen(clsName, mg, cp, factory, il, argsMap, argsStartPos, funcRefsMap);
		arg2.bytecodeGen(clsName, mg, cp, factory, il, argsMap, argsStartPos, funcRefsMap);
		return il.append(InstructionConstants.DADD);
	}

}
