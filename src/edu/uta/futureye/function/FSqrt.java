package edu.uta.futureye.function;

import static com.sun.org.apache.bcel.internal.generic.InstructionConstants.DALOAD;

import java.util.Map;

import com.sun.org.apache.bcel.internal.Constants;
import com.sun.org.apache.bcel.internal.generic.ALOAD;
import com.sun.org.apache.bcel.internal.generic.ConstantPoolGen;
import com.sun.org.apache.bcel.internal.generic.InstructionFactory;
import com.sun.org.apache.bcel.internal.generic.InstructionHandle;
import com.sun.org.apache.bcel.internal.generic.InstructionList;
import com.sun.org.apache.bcel.internal.generic.MethodGen;
import com.sun.org.apache.bcel.internal.generic.PUSH;
import com.sun.org.apache.bcel.internal.generic.Type;

import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.util.Constant;

public class FSqrt extends AbstractSimpleMathFunc {

	public FSqrt() {
		super("sqrt", Constant.x);
	}
	
	public FSqrt(String varName) {
		super("sqrt", varName);
	}

	@Override
	public double apply(double... args) {
		return Math.sqrt(args[argIdx]);
	}

	@Override
	public double apply(Variable v) {
		return Math.sqrt(v.get(varName));
	}
	@Override
	public MathFunc diff(String varName) {
		return new FPow(0.5).diff(varName).setArgIdx(this.getArgIdxMap());
	}
	
	@Override
	public InstructionHandle bytecodeGen(String clsName, MethodGen mg, 
			ConstantPoolGen cp, InstructionFactory factory, 
			InstructionList il, Map<String, Integer> argsMap, 
			int argsStartPos, Map<MathFunc, Integer> funcRefsMap) {
		il.append(new ALOAD(argsStartPos));
		il.append(new PUSH(cp, argsMap.get(varName)));
		il.append(DALOAD);
		return  il.append(factory.createInvoke("java.lang.Math", "sqrt",
				Type.DOUBLE, new Type[] { Type.DOUBLE }, 
				Constants.INVOKESTATIC));
	}
	
	@Override
	public String getExpr() {
		return fName + "(" + varName + ")";
	}
	
	@Override
	public String toString() {
		return fName + "(" + varName + ")";
	}
	
}
