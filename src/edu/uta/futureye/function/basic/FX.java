package edu.uta.futureye.function.basic;

import static org.apache.bcel.generic.InstructionConstants.DALOAD;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.bcel.generic.ALOAD;
import org.apache.bcel.generic.ConstantPoolGen;
import org.apache.bcel.generic.InstructionFactory;
import org.apache.bcel.generic.InstructionHandle;
import org.apache.bcel.generic.InstructionList;
import org.apache.bcel.generic.MethodGen;
import org.apache.bcel.generic.PUSH;
import org.objectweb.asm.MethodVisitor;

import com.sun.xml.internal.ws.org.objectweb.asm.Opcodes;

import edu.uta.futureye.core.Node;
import edu.uta.futureye.function.FMath;
import edu.uta.futureye.function.MathFuncBasic;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.util.Constant;

/**
 * Identity function
 * f(x) = x
 * 
 */
public class FX extends MathFuncBasic {
	/**
	 * Predefined instances of FX
	 */
	public final static FX x = new FX(Constant.x); 
	public final static FX y = new FX(Constant.y); 
	public final static FX z = new FX(Constant.z); 
	
	public final static FX r = new FX(Constant.r); 
	public final static FX s = new FX(Constant.s); 
	public final static FX t = new FX(Constant.t); 
	
	String varName;
	int argIdx; //the index of varName in args of method apply(double ..args)
	
	/**
	 * Identity function: f(varName) = varName
	 */
	public FX(String varName) {
		this.varName = varName;
		this.argIdx = 0;
	}

	/**
	 * This is used only in non-compiled version
	 * argIdx is introduced for this function with argument 'double ...args'
	 * argIdx is not used in compiled version instead it uses 'argsMap'
	 * see bytecodeGen(...,argsMap,...)
	 * 
	 */
	@Override
	public double apply(double... args) {
		return args[argIdx];
	}
	
	/**
	 * This version doesn't have the problem in finding the index
	 * in the version 'apply(double... args)'
	 */
	@Deprecated
	@Override
	public double apply(Variable v) {
		return v.get(varName);
	}

	@Override
	public MathFunc diff(String varName) {
		if(this.varName.equals(varName))
			return FMath.C1;
		else
			return FMath.C0;
	}
	
	@Override
	public String getExpr() {
		return varName;
	}
	
	@Override
	public String toString() {
		return varName;
	}

	@Override
	public InstructionHandle bytecodeGen(String clsName, MethodGen mg,
			ConstantPoolGen cp, InstructionFactory factory,
			InstructionList il, Map<String, Integer> argsMap, int argsStartPos, 
			Map<MathFunc, Integer> funcRefsMap) {
		il.append(new ALOAD(argsStartPos));
		il.append(new PUSH(cp, argsMap.get(varName)));
		return il.append(DALOAD);
	}

	@Override
	public void bytecodeGen(MethodVisitor mv, Map<String, Integer> argsMap,
			int argsStartPos, Map<MathFunc, Integer> funcRefsMap, String clsName) {
		//No need to compile to static field for performance consideration
//		if (this.compileToStaticField && !this.isCompiledToStaticFiled) {
//			mv.visitFieldInsn(Opcodes.GETSTATIC, this.genClassName, this.staticFieldName, "D");
//			this.isCompiledToStaticFiled = true;
//		} else {
			mv.visitIntInsn(Opcodes.ALOAD, argsStartPos);
			mv.visitLdcInsn(argsMap.get(varName));
			mv.visitInsn(Opcodes.DALOAD);
//		}
	}
	
	@Override
	public MathFunc setName(String name) {
		this.varName = name;
		return this;
	}

	@Override
	public MathFunc setVarNames(List<String> varNames) {
		this.varName = varNames.get(0);
		return this;
	}
	
	@Override
	public List<String> getVarNames() {
		List<String> ret = new ArrayList<String>();
		ret.add(varName);
		return ret;
	}

	public String getVarName() {
		return this.varName;
	}
	
	public MathFunc setVarName(String varName) {
		this.varName = varName;
		return this;
	}

	@Override
	public MathFunc setArgIdx(Map<String, Integer> argsMap) {
		this.argIdx = argsMap.get(varName);
		return this;
	}
	
	@Override
	public Map<String, Integer> getArgIdxMap() {
		Map<String, Integer> ret = new HashMap<String, Integer>();
		ret.put(varName, argIdx);
		return ret;
	}
	
	@Override
	public boolean isConstant() {
		return false;
	}
	
	@Override
	public boolean isInteger() {
		return false;
	}
	
	@Override
	public boolean isZero() {
		return false;
	}
	
	@Override
	public boolean isReal() {
		return false;
	}

}
