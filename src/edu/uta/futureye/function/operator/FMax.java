package edu.uta.futureye.function.operator;

import java.util.Map;

import org.apache.bcel.Constants;
import org.apache.bcel.generic.ConstantPoolGen;
import org.apache.bcel.generic.InstructionFactory;
import org.apache.bcel.generic.InstructionHandle;
import org.apache.bcel.generic.InstructionList;
import org.apache.bcel.generic.MethodGen;
import org.apache.bcel.generic.Type;

import edu.uta.futureye.function.intf.MathFunc;

public class FMax extends FBinaryOp {
	/**
	 * Construct function: max(f(x), g(y))
	 * 
	 * @param f
	 * @param g
	 */
	public FMax(MathFunc f, MathFunc g) {
		super(f, g);
	}
	
	@Override
	public double apply(double... args) {
		return Math.max(arg1.apply(args), arg2.apply(args));
	}
	
	@Override
	public InstructionHandle bytecodeGen(String clsName, MethodGen mg,
			ConstantPoolGen cp, InstructionFactory factory,
			InstructionList il, Map<String, Integer> argsMap, int argsStartPos, 
			Map<MathFunc, Integer> funcRefsMap) {
		arg1.bytecodeGen(clsName, mg, cp, factory, il, argsMap, argsStartPos, funcRefsMap);
		arg2.bytecodeGen(clsName, mg, cp, factory, il, argsMap, argsStartPos, funcRefsMap);
		return  il.append(factory.createInvoke("java.lang.Math", "max",
				Type.DOUBLE, 
				new Type[] { Type.DOUBLE, Type.DOUBLE },
		Constants.INVOKESTATIC));
	}
	
	@Override
	public String getExpr() {
		return "max("+arg1.getExpr()+","+arg2.getExpr()+")";
	}
}
