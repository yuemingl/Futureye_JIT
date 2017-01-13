package edu.uta.futureye.function.operator;

import java.util.Map;

import org.apache.bcel.Constants;
import org.apache.bcel.generic.ConstantPoolGen;
import org.apache.bcel.generic.InstructionFactory;
import org.apache.bcel.generic.InstructionHandle;
import org.apache.bcel.generic.InstructionList;
import org.apache.bcel.generic.MethodGen;
import org.apache.bcel.generic.Type;
import org.objectweb.asm.MethodVisitor;

import com.sun.xml.internal.ws.org.objectweb.asm.Opcodes;

import edu.uta.futureye.function.FMath;
import edu.uta.futureye.function.intf.MathFunc;

public class FLog10 extends FUniaryOp {
	/**
	 * Construct function: log10(g(x))
	 * 
	 * @param base
	 * @param g
	 */
	public FLog10(MathFunc g) {
		super(g);
	}
	
	@Override
	public double apply(double... args) {
		return Math.log10(arg.apply(args));
	}

	@Override
	public MathFunc diff(String varName) {
		MathFunc ret = FMath.C1.D(this.M(Math.log(10))).M(arg.diff(varName));
		return ret.setArgIdx(this.getArgIdxMap());
	}
	
	@Override
	public InstructionHandle bytecodeGen(String clsName, MethodGen mg,
			ConstantPoolGen cp, InstructionFactory factory,
			InstructionList il, Map<String, Integer> argsMap, int argsStartPos, 
			Map<MathFunc, Integer> funcRefsMap) {
		arg.bytecodeGen(clsName, mg, cp, factory, il, argsMap, argsStartPos, funcRefsMap);
		return  il.append(factory.createInvoke("java.lang.Math", "log10",
				Type.DOUBLE, 
				new Type[] { Type.DOUBLE },
		Constants.INVOKESTATIC));
	}
	
	@Override
	public void bytecodeGen(MethodVisitor mv, Map<String, Integer> argsMap,
			int argsStartPos, Map<MathFunc, Integer> funcRefsMap, String clsName) {
		if (this.compileToStaticField && this.isCompiledToStaticFiled) {
			mv.visitFieldInsn(Opcodes.GETSTATIC, this.genClassName, this.staticFieldName, "D");
		} else {
			arg.bytecodeGen(mv, argsMap, argsStartPos, funcRefsMap, clsName);
			mv.visitMethodInsn(Opcodes.INVOKESTATIC, "java/lang/Math", "log10", "(D)D", false);
		}
	}
	
	@Override
	public String getExpr() {
		return "log10("+arg.getExpr()+")";
	}
}
