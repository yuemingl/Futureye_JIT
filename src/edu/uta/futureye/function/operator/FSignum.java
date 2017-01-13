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

import edu.uta.futureye.function.intf.MathFunc;

public class FSignum extends FUniaryOp {
	/**
	 * Construct function : signum(g(x))
	 * 
	 * @param g
	 */
	public FSignum(MathFunc g) {
		super(g);
	}

	@Override
	public double apply(double... args) {
		return Math.signum(arg.apply(args));
	}
	
	@Override
	public InstructionHandle bytecodeGen(String clsName, MethodGen mg,
			ConstantPoolGen cp, InstructionFactory factory,
			InstructionList il, Map<String, Integer> argsMap, int argsStartPos, 
			Map<MathFunc, Integer> funcRefsMap) {
		arg.bytecodeGen(clsName, mg, cp, factory, il, argsMap, argsStartPos, funcRefsMap);
		return  il.append(factory.createInvoke("java.lang.Math", "signum",
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
			mv.visitMethodInsn(Opcodes.INVOKESTATIC, "java/lang/Math", "signum", "(D)D", false);
		}
	}

	@Override
	public String getExpr() {
		return "signum("+arg.getExpr()+")";
	}
}
