package edu.uta.futureye.function.operator;

import java.util.Map;

import com.sun.org.apache.bcel.internal.Constants;
import com.sun.org.apache.bcel.internal.generic.ConstantPoolGen;
import com.sun.org.apache.bcel.internal.generic.InstructionFactory;
import com.sun.org.apache.bcel.internal.generic.InstructionHandle;
import com.sun.org.apache.bcel.internal.generic.InstructionList;
import com.sun.org.apache.bcel.internal.generic.MethodGen;
import com.sun.org.apache.bcel.internal.generic.PUSH;
import com.sun.org.apache.bcel.internal.generic.Type;

import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.intf.MathFunc;

public class FPow2 extends FBinaryOp {
	/**
	 * 
	 * Construct function: pow(base(x), exp(y))
	 * 
	 * @param base
	 * @param exp
	 */
	public FPow2(MathFunc base, MathFunc exp) {
		super(base, exp);
	}
	
	public FPow2(MathFunc base, int exp) {
		super(base, FC.c(exp));
	}
	
	@Override
	public double apply(double... args) {
		throw new UnsupportedOperationException();
	}
	
	@Override
	public InstructionHandle bytecodeGen(String clsName, MethodGen mg,
			ConstantPoolGen cp, InstructionFactory factory,
			InstructionList il, Map<String, Integer> argsMap, int argsStartPos, 
			Map<MathFunc, Integer> funcRefsMap) {
		arg1.bytecodeGen(clsName, mg, cp, factory, il, argsMap, argsStartPos, funcRefsMap);
		if(arg2 instanceof FC) {
			il.append(new PUSH(cp, (int)arg2.apply()));
			return  il.append(factory.createInvoke("edu.uta.futureye.function.operator.FPow2", "powi",
					Type.DOUBLE, 
					new Type[] { Type.DOUBLE, Type.INT },
			Constants.INVOKESTATIC));
		} else {
			arg2.bytecodeGen(clsName, mg, cp, factory, il, argsMap, argsStartPos, funcRefsMap);
			return  il.append(factory.createInvoke("java.lang.Math", "pow",
					Type.DOUBLE, 
					new Type[] { Type.DOUBLE, Type.DOUBLE },
			Constants.INVOKESTATIC));
		}
	}
	
	@Override
	public String getExpr() {
		return "pow("+arg1.getExpr()+","+arg2.getExpr()+")";
	}
	
	public static double powi(double base, int exp) {
		if(exp == 0) return 1.0;
		else if(exp < 0) return 1.0/powi(base, -exp);
		else if(exp == 1) return base;
		else { //exp >= 2
			double rlt = 1.0;
			double tmp = base;
			if((exp & 0x1) > 0) rlt = base;
			int mask = exp>>>1;
			while(mask > 0) {
				tmp *= tmp;
				if((mask & 0x1) > 0) {
					rlt *= tmp;
				}
				mask >>>= 1;
			}
			return rlt;
		}
	}
}
