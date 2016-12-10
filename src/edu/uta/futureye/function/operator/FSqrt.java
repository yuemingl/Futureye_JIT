package edu.uta.futureye.function.operator;

import java.util.Map;

import org.apache.bcel.Constants;
import org.apache.bcel.generic.ConstantPoolGen;
import org.apache.bcel.generic.InstructionFactory;
import org.apache.bcel.generic.InstructionHandle;
import org.apache.bcel.generic.InstructionList;
import org.apache.bcel.generic.MethodGen;
import org.apache.bcel.generic.Type;

import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.intf.MathFunc;

public class FSqrt extends FUniaryOp {
	/**
	 * Construct function : sqrt(g(x))
	 * 
	 * @param g
	 */
	public FSqrt(MathFunc g) {
		super(g);
	}

	@Override
	public double apply(double... args) {
		return Math.sqrt(arg.apply(args));
	}
	
	@Override
	public MathFunc diff(String varName) {
		MathFunc ret = FC.c(0.5).M(new FPow(arg,FC.c(-0.5))).M(arg.diff(varName));
		return ret.setArgIdx(this.getArgIdxMap());
	}
	
	@Override
	public InstructionHandle bytecodeGen(String clsName, MethodGen mg,
			ConstantPoolGen cp, InstructionFactory factory,
			InstructionList il, Map<String, Integer> argsMap, int argsStartPos, 
			Map<MathFunc, Integer> funcRefsMap) {
		arg.bytecodeGen(clsName, mg, cp, factory, il, argsMap, argsStartPos, funcRefsMap);
		return  il.append(factory.createInvoke("java.lang.Math", "sqrt",
				Type.DOUBLE, 
				new Type[] { Type.DOUBLE },
		Constants.INVOKESTATIC));
	}
	
	@Override
	public String getExpr() {
		return "sqrt("+arg.getExpr()+")";
	}
}
