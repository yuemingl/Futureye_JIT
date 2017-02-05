package edu.uta.futureye.function.operator;

import java.util.Map;

import org.apache.bcel.generic.ConstantPoolGen;
import org.apache.bcel.generic.InstructionConstants;
import org.apache.bcel.generic.InstructionFactory;
import org.apache.bcel.generic.InstructionHandle;
import org.apache.bcel.generic.InstructionList;
import org.apache.bcel.generic.MethodGen;
import org.objectweb.asm.MethodVisitor;

import com.sun.xml.internal.ws.org.objectweb.asm.Opcodes;

import edu.uta.futureye.function.VariableArray;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.lib.assembler.AssembleParam;
import edu.uta.futureye.util.Constant;

public class FSub extends FBinaryOp {
	public FSub(MathFunc left, MathFunc right) {
		super(left, right);
	}

	@Override
	public double apply(double... args) {
		return arg1.apply(args) - arg2.apply(args);
	}
	
	@Override
	public double apply(AssembleParam ap, double... args) {
		return arg1.apply(ap,args) - arg2.apply(ap, args);
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
	public String getExpr() {
		StringBuilder sb = new StringBuilder();
		if(! (arg1.isConstant() && Math.abs(arg1.apply()) < Constant.eps) ) {
			sb.append(arg1.getExpr());
		}
		sb.append(" - ");
		if(arg2.getOpOrder() >= OP_ORDER3)
			sb.append("(").append(arg2.getExpr()).append(")");
		else
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
		return il.append(InstructionConstants.DSUB);
	}

	@Override
	public void bytecodeGen(MethodVisitor mv, Map<String, Integer> argsMap,
			int argsStartPos, Map<MathFunc, Integer> funcRefsMap, String clsName) {
		if (this.compileToStaticField && this.isCompiledToStaticFiled) {
			mv.visitFieldInsn(Opcodes.GETSTATIC, this.genClassName, this.staticFieldName, "D");
		} else {
			arg1.bytecodeGen(mv, argsMap, argsStartPos, funcRefsMap, clsName);
			arg2.bytecodeGen(mv, argsMap, argsStartPos, funcRefsMap, clsName);
			mv.visitInsn(Opcodes.DSUB);
		}
	}

}
