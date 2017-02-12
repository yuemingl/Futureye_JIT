package edu.uta.futureye.function.basic;

import java.util.Map;

import org.apache.bcel.generic.ALOAD;
import org.apache.bcel.generic.ConstantPoolGen;
import org.apache.bcel.generic.DALOAD;
import org.apache.bcel.generic.InstructionFactory;
import org.apache.bcel.generic.InstructionHandle;
import org.apache.bcel.generic.InstructionList;
import org.apache.bcel.generic.MethodGen;
import org.apache.bcel.generic.PUSH;
import org.objectweb.asm.MethodVisitor;

import com.sun.xml.internal.ws.org.objectweb.asm.Opcodes;

import edu.uta.futureye.function.FMath;
import edu.uta.futureye.function.SingleVarFunc;
import edu.uta.futureye.function.intf.MathFunc;


public class TriAreaCoordS extends SingleVarFunc {
	MathFunc jac;
	MathFunc x1;
	MathFunc x2;
	MathFunc x3;
	MathFunc y1;
	MathFunc y2;
	MathFunc y3;

	public TriAreaCoordS() {
		super("s", "s");
	}

	@Override
	public double apply(double... args) {
		return args[this.argIdx];
	}

	public void setJacXY(MathFunc jac, 
			MathFunc x1, MathFunc x2, MathFunc x3, 
			MathFunc y1, MathFunc y2, MathFunc y3 ) {
		this.jac = jac;
		this.x1 = x1;
		this.x2 = x2;
		this.x3 = x3;
		this.y1 = y1;
		this.y2 = y2;
		this.y3 = y3;
	}
	
//	s_x = (y3-y1)/jac;
//	s_y = (x1-x3)/jac;
	@Override
	public MathFunc diff(String varName) {
		if(varName.equals("s"))
			return FMath.C1;
		if(varName.equals("x"))
			return (y3-y1)/jac;
		else if(varName.equals("y"))
			return (x1-x3)/jac;
		else
			return FMath.C0;
	}
	
	public String getExpr() {
		return this.varName;
	}
	
	@Override
	public void bytecodeGen(MethodVisitor mv, Map<String, Integer> argsMap,
			int argsStartPos, Map<MathFunc, Integer> funcRefsMap,
			String clsName) {
		mv.visitIntInsn(Opcodes.ALOAD, argsStartPos);
		mv.visitLdcInsn(argsMap.get(varName));
		mv.visitInsn(Opcodes.DALOAD);
	}
	
	@Override
	public InstructionHandle bytecodeGen(String clsName, MethodGen mg, 
			ConstantPoolGen cp, InstructionFactory factory, 
			InstructionList il, Map<String, Integer> argsMap, 
			int argsStartPos, Map<MathFunc, Integer> funcRefsMap) {
		il.append(new ALOAD(argsStartPos));
		il.append(new PUSH(cp, argsMap.get(this.getName())));
		return il.append(new DALOAD());
	}
}