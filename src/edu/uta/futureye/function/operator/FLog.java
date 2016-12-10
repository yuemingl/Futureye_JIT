package edu.uta.futureye.function.operator;

import java.util.Map;

import org.apache.bcel.Constants;
import org.apache.bcel.generic.ConstantPoolGen;
import org.apache.bcel.generic.InstructionFactory;
import org.apache.bcel.generic.InstructionHandle;
import org.apache.bcel.generic.InstructionList;
import org.apache.bcel.generic.MethodGen;
import org.apache.bcel.generic.Type;

import edu.uta.futureye.function.FMath;
import edu.uta.futureye.function.intf.MathFunc;

/**
 * The natural logarithm (base e) of an instance of MathFunc
 *
 */
public class FLog extends FUniaryOp {
	/**
	 * Construct function: log(g(x))
	 * 
	 * @param base
	 * @param g
	 */
	public FLog(MathFunc g) {
		super(g);
	}
	
	@Override
	public double apply(double... args) {
		return Math.log(arg.apply(args));
	}

	@Override
	public MathFunc diff(String varName) {
		MathFunc ret = FMath.C1.D(arg).M(arg.diff(varName));
		return ret.setArgIdx(this.getArgIdxMap());
	}
	
	@Override
	public InstructionHandle bytecodeGen(String clsName, MethodGen mg,
			ConstantPoolGen cp, InstructionFactory factory,
			InstructionList il, Map<String, Integer> argsMap, int argsStartPos, 
			Map<MathFunc, Integer> funcRefsMap) {
		arg.bytecodeGen(clsName, mg, cp, factory, il, argsMap, argsStartPos, funcRefsMap);
		return  il.append(factory.createInvoke("java.lang.Math", "log",
				Type.DOUBLE, 
				new Type[] { Type.DOUBLE },
		Constants.INVOKESTATIC));
	}
	
	@Override
	public String getExpr() {
		return "log("+arg.getExpr()+")";
	}
}
