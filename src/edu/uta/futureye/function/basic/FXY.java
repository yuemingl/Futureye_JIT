package edu.uta.futureye.function.basic;

import java.util.List;
import java.util.Map;
import java.util.UUID;

import org.apache.bcel.generic.ConstantPoolGen;
import org.apache.bcel.generic.InstructionFactory;
import org.apache.bcel.generic.InstructionHandle;
import org.apache.bcel.generic.InstructionList;
import org.apache.bcel.generic.MethodGen;
import org.objectweb.asm.MethodVisitor;

import com.sun.xml.internal.ws.org.objectweb.asm.Opcodes;

import edu.uta.futureye.function.FMath;
import edu.uta.futureye.function.MultiVarFunc;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.util.Constant;

/**
 * f(x,y) = c1*x + c2*y + c3
 * 
 */
public class FXY extends MultiVarFunc{
	protected double c1, c2, c3;

	public FXY(double c1, double c2, double c3) {
		this.varNames = new String[2];
		this.varNames[0] = Constant.x;
		this.varNames[1] = Constant.y;
		this.c1 = c1;
		this.c2 = c2;
		this.c3 = c3;
	}	
	
	public FXY(String funcName, List<String> varNames,double c1,double c2,double c3) {
		super(funcName, varNames);
		this.c1 = c1;
		this.c2 = c2;
		this.c3 = c3;
	}

	@Override
	public double apply(double... args) {
		return c1*args[this.argIdx[0]] + c2*args[this.argIdx[1]] + c3;
	}
	
	@Override
	public MathFunc diff(String varName) {
		if(varNames[0].equals(varName))
			return new FC(c1);
		else if(varNames[1].equals(varName)) {
			return new FC(c2);
		}
		return FMath.C0;
	}

	@Override
	public InstructionHandle bytecodeGen(String clsName, MethodGen mg,
			ConstantPoolGen cp, InstructionFactory factory,
			InstructionList il, Map<String, Integer> argsMap, int argsStartPos, 
			Map<MathFunc, Integer> funcRefsMap) {
		throw new UnsupportedOperationException();
	}

	@Override
	public void bytecodeGen(MethodVisitor mv, Map<String, Integer> argsMap,
			int argsStartPos, Map<MathFunc, Integer> funcRefsMap, String clsName) {
		mv.visitIntInsn(Opcodes.ALOAD, argsStartPos);
		mv.visitLdcInsn(argsMap.get(this.varNames[0]));
		mv.visitInsn(Opcodes.DALOAD);
		mv.visitLdcInsn(this.c1);
		mv.visitInsn(Opcodes.DMUL);
		mv.visitIntInsn(Opcodes.ALOAD, argsStartPos);
		mv.visitLdcInsn(argsMap.get(this.varNames[1]));
		mv.visitInsn(Opcodes.DALOAD);
		mv.visitLdcInsn(this.c2);
		mv.visitInsn(Opcodes.DMUL);
		mv.visitInsn(Opcodes.DADD);
		mv.visitLdcInsn(this.c3);
		mv.visitInsn(Opcodes.DADD);
	}
	
	public String toString() {
		//fixme c2=0
		String s1 = "";
		if(Math.abs(c1) > Constant.eps)
			s1 = c1 + "*" + this.varNames[0];
		String s2 = "";
		if(Math.abs(c2) > Constant.eps)
			s2 = c2 + "*" + this.varNames[1];
		String p1 = (s1.isEmpty()||s2.isEmpty())?"":" + ";
		String s3 = "";
		if(Math.abs(c3) > Constant.eps)
			s3 = c3;
		String p2 = (s2.isEmpty()||s3.isEmpty())?"":" + ";
		return s1+p1+s2+p2+s3;
	}
	
	public String getExpr() {
		return toString();
	}
	
	@Override
	public String getName() {
		return this.fName;
	}
}