/**
 * Copyright (c) 2010, nkliuyueming@gmail.com. All rights reserved.
 * 
 * 
 */
package edu.uta.futureye.function;

import static org.apache.bcel.Constants.ACC_PUBLIC;

import java.io.File;
import java.io.FileOutputStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.bcel.Constants;
import org.apache.bcel.generic.AALOAD;
import org.apache.bcel.generic.ALOAD;
import org.apache.bcel.generic.ArrayType;
import org.apache.bcel.generic.ClassGen;
import org.apache.bcel.generic.ConstantPoolGen;
import org.apache.bcel.generic.FieldGen;
import org.apache.bcel.generic.GETFIELD;
import org.apache.bcel.generic.InstructionFactory;
import org.apache.bcel.generic.InstructionHandle;
import org.apache.bcel.generic.InstructionList;
import org.apache.bcel.generic.MethodGen;
import org.apache.bcel.generic.PUSH;
import org.apache.bcel.generic.Type;
import org.objectweb.asm.FieldVisitor;
import org.objectweb.asm.Label;
import org.objectweb.asm.MethodVisitor;
import org.objectweb.asm.Opcodes;

import edu.uta.futureye.bytecode.CompiledFunc;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.basic.FComposite;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.function.operator.FAdd;
import edu.uta.futureye.function.operator.FDiv;
import edu.uta.futureye.function.operator.FMul;
import edu.uta.futureye.function.operator.FSub;
import edu.uta.futureye.lib.assembler.AssembleParam;
import edu.uta.futureye.util.BytecodeConst;
import edu.uta.futureye.util.BytecodeUtils;
import edu.uta.futureye.util.ClassGenerator;
import edu.uta.futureye.util.Constant;
import edu.uta.futureye.util.FuncClassLoader;

public abstract class MathFuncBase implements MathFunc, Cloneable { 
	
	@Override
	public abstract double apply(double ...args);
	
	@Override
	public double apply(AssembleParam ap, double... args) {
		return apply(args);
	}

	@Deprecated
	@Override
	public double apply(Variable v) {
		Node n = new Node(v.getIndex());
		return apply(new AssembleParam(v.getElement(),-1,-1), v.getVarValues());
	}
	
	@Deprecated
	@Override
	public double apply(Variable v, Map<Object,Object> cache) {
		//Ignore cache
		return apply(v);
	}
	
	@Deprecated
	@Override
	public double[] applyAll(VariableArray v, Map<Object,Object> cache) {
		throw new UnsupportedOperationException();
	}
	
	@Override
	public MathFunc diff(String varName) {
		throw new UnsupportedOperationException();
	}
	
	@Override
	public MathFunc compose(final Map<String,MathFunc> fInners) {
		boolean find = false;
		for(String key : fInners.keySet()) {
			if(getVarNames().contains(key)) find = true;
		}
		if(!find) 
			return this; //No compose
		else {
			MathFunc ret = new FComposite(this, fInners);
//this is not have to be a composite function. so, this is useless
//			if(this instanceof FComposite && this.isOuterVarActive())
//				ret.setOuterVarActive();
//			else
//				ret.setInnerVarActive();
			return ret;
		}
	}
	
	////////////////////////Basic Math Operations/////////////////////////////
	
	@Override
	public MathFunc A(MathFunc g) {
		final MathFunc f1 = this;
		final MathFunc f2 = g;
		if(f1.isConstant() && f2.isConstant()) {
			return new FC(f1.apply() + f2.apply());
		} else if(f1.isConstant() && Math.abs(f1.apply()) < Constant.eps) {
			return f2;
		} else if(f2.isConstant() && Math.abs(f2.apply()) < Constant.eps) {
			return f1;
		} else {
			return new FAdd(f1, f2);
		}
	}
	@Override
	public MathFunc A(double g) {
		return A(FC.c(g));
	}
	
	@Override
	public MathFunc S(MathFunc g) {
		final MathFunc f1 = this;
		final MathFunc f2 = g;
		if(f1.isConstant() && f2.isConstant()) {
			return new FC(f1.apply() - f2.apply());
		} else if(f2.isConstant() && Math.abs(f2.apply()) < Constant.eps) {
			return f1;
		} else {
			return new FSub(f1, f2);
		}
	}
	@Override
	public MathFunc S(double g) {
		return S(FC.c(g));
	}	
	
	@Override
	public MathFunc M(MathFunc f) {
		final MathFunc f1 = this;
		final MathFunc f2 = f;
		if(f1.isConstant() && f2.isConstant()) {
			return new FC(f1.apply() * f2.apply());
		} else if( (f1.isConstant() && Math.abs(f1.apply()) < Constant.eps) ||
				f2.isConstant() && Math.abs(f2.apply()) < Constant.eps)
			return FMath.C0;
		else if(f1.isConstant() && Math.abs(f1.apply()-1.0) < Constant.eps)
			return f2;
		else if(f2.isConstant() && Math.abs(f2.apply()-1.0) < Constant.eps)
			return f1;
		else
			return new FMul(f1, f2);
	}
	@Override
	public MathFunc M(double g) {
		return M(FC.c(g));
	}	
	
	@Override
	public MathFunc D(MathFunc f) {
		final MathFunc f1 = this;
		final MathFunc f2 = f;
		if(f1.isConstant() && f2.isConstant()) {
			return new FC(f1.apply() / f2.apply());
		} else if(f1.isConstant() && Double.compare(f1.apply(),0.0)==0) {
			//Math.abs(f1.value())<Constant.eps will not work properly
			return FMath.C0;
		} else if(f2.isConstant() && Double.compare(f2.apply(),0.0)==0) {
			return FC.c(Double.POSITIVE_INFINITY);
		}  else if(f2.isConstant() && Math.abs(f2.apply()-1.0) < Constant.eps) {
			return f1;
		} else {
			return new FDiv(f1, f2);
		}
	}
	@Override
	public MathFunc D(double g) {
		return D(FC.c(g));
	}
	
	///////////////////////////Compilation///////////////////////////////
	
	/**
	 * This implementation calls the 'apply(...)' method in referenced function by default 
	 * in the case of that a generated sub-class does not override this method (bytecodeGen)
	 */
	@Override
	public InstructionHandle bytecodeGen(String clsName, MethodGen mg, 
			ConstantPoolGen cp, InstructionFactory factory, 
			InstructionList il, Map<String, Integer> argsMap, 
			int argsStartPos, Map<MathFunc, Integer> funcRefsMap) {
//		throw new UnsupportedOperationException();
		FieldGen fg = new FieldGen(ACC_PUBLIC, new ArrayType(Type.getType(MathFunc.class), 1), "funcRefs", cp);
		//System.out.println(fg.getSignature());
		int idxFuncRefs = cp.addFieldref("edu.uta.futureye.bytecode.CompiledFunc", "funcRefs", fg.getSignature());

		il.append(new ALOAD(0));
		il.append(new GETFIELD(idxFuncRefs));
		il.append(new PUSH(cp, funcRefsMap.get(this)));
		il.append(new AALOAD());
		il.append(new ALOAD(BytecodeConst.assembleParamIdx+1));
		il.append(new ALOAD(BytecodeConst.argIdx+1));
		return  il.append(factory.createInvoke("edu.uta.futureye.function.intf.MathFunc", "apply",
				Type.DOUBLE, new Type[] {
					Type.getType(AssembleParam.class),
					new ArrayType(Type.DOUBLE, 1)
				}, 
				Constants.INVOKEINTERFACE));
	}

	@Override
	public void bytecodeGen(MethodVisitor mv, Map<String, Integer> argsMap,
			int argsStartPos, Map<MathFunc, Integer> funcRefsMap,
			String clsName) {
		mv.visitVarInsn(org.objectweb.asm.Opcodes.ALOAD, 0);
		mv.visitFieldInsn(Opcodes.GETFIELD, ClassGenerator.getASMName(CompiledFunc.class), "funcRefs", 
				org.objectweb.asm.Type.getType(MathFunc[].class).getDescriptor());
		mv.visitLdcInsn(funcRefsMap.get(this));
		mv.visitInsn(org.objectweb.asm.Opcodes.AALOAD);
		mv.visitVarInsn(org.objectweb.asm.Opcodes.ALOAD, BytecodeConst.assembleParamIdx+1);
		mv.visitVarInsn(org.objectweb.asm.Opcodes.ALOAD, BytecodeConst.argIdx+1);
		mv.visitMethodInsn(org.objectweb.asm.Opcodes.INVOKEINTERFACE, 
				ClassGenerator.getASMName(MathFunc.class), "apply", "("+
						org.objectweb.asm.Type.getType(AssembleParam.class).getDescriptor()+
						org.objectweb.asm.Type.getType(double[].class).getDescriptor()+
						")D", true);
	}
	
	@Override
	public CompiledFunc compile(String ...varNames) {
		String clsName = getName();
		if(clsName == null || clsName.length() == 0)
			clsName = this.getClass().getSimpleName();
		clsName = clsName + java.util.UUID.randomUUID().toString().replaceAll("-", "");
		
		FuncClassLoader<CompiledFunc> fcl = new FuncClassLoader<CompiledFunc>(CompiledFunc.class.getClassLoader());
		ClassGen genClass = BytecodeUtils.genClass(this, varNames, clsName, true, false);
		CompiledFunc func = fcl.newInstance(genClass);
		
		// Set funcRefs field in CompiledFunc
		List<MathFunc> list = new ArrayList<MathFunc>();
		BytecodeUtils.postOrder(this, list);
		func.setFuncRefs(list.toArray(new MathFunc[0]));
		
		return func;
	}

	@Override
	public CompiledFunc compileWithASM(String ...varNames) {

		boolean writeFile = true;
		genClassName = getName();
		if (genClassName == null || genClassName.length() == 0)
			genClassName = this.getClass().getSimpleName();
		genClassName = genClassName
				+ java.util.UUID.randomUUID().toString().replaceAll("-", "");
		try {
			FuncClassLoader<CompiledFunc> mcl = FuncClassLoader.getInstance(ClassGenerator.class.getClassLoader());
			ClassGenerator cgen = new ClassGenerator(genClassName);

			cgen.startClass(ClassGenerator.getASMName(CompiledFunc.class), null);

			// Define method:
			// double apply(Element e, Node n, double ...args);
			MethodVisitor mv = null;
			String methodName = "apply";

			// Generate the function for the root expression
			Label startMatchesLabel = new Label();
			Label endMatchesLabel = new Label();
			org.objectweb.asm.Type retType = org.objectweb.asm.Type.getType(double.class);
			org.objectweb.asm.Type param1 = org.objectweb.asm.Type.getType(AssembleParam.class);
			org.objectweb.asm.Type param2 = org.objectweb.asm.Type.getType(double[].class);
			mv = cgen.startMethod(Opcodes.ACC_PUBLIC, methodName,
					org.objectweb.asm.Type.getMethodDescriptor(retType, param1, param2));
			cgen.startCode(mv, startMatchesLabel);

			HashMap<String, Integer> argsMap = new HashMap<String, Integer>();
			if(varNames == null || varNames.length == 0) {
				List<String> args = this.getVarNames();
				for(int i=0; i<args.size(); i++) {
					argsMap.put(args[i], i);
				}
				//System.out.println("JIT compileWithASM(): "+this.getName()+"("+argsMap+")"+" = "+this.getExpr());
			} else {
				StringBuilder sb = new StringBuilder();
				sb.append("(");
				for(int i=0; i<varNames.length; i++) {
					argsMap.put(varNames[i], i);
					sb.append(varNames[i]).append(",");
				}
				sb.delete(sb.length()-1, sb.length());
				sb.append(")");
				//func.setArgIdx(argsMap); //No need, this is for user defined apply() method not for compile()
				String expr = this.getExpr();
				if(expr.length()>50)
					expr = expr.substring(0,50);
				//System.out.println("JIT compileWithASM: "+this.getName()+sb.toString()+" = "+expr);
			}

			Map<MathFunc, Integer> refsMap = BytecodeUtils.getFuncRefsMap(this);
			
			if (this.compileToStaticField) {
				this.bytecodeGen(mv, argsMap, 2, refsMap, genClassName); //2 for args: double apply(Element e, Node n, double ...args);
				staticFieldName = "var_" + this.getName();
				FieldVisitor fv = cgen.getClassWriter().visitField(
						Opcodes.ACC_PUBLIC + Opcodes.ACC_STATIC, staticFieldName, "D", "D", 0.0);
				fv.visitEnd();
				mv.visitInsn(Opcodes.DUP2);
				mv.visitFieldInsn(Opcodes.PUTSTATIC, cgen.getClassName(),
						staticFieldName, "D");
				
				//no need to do this here:
				//this.compileToStaticField = false; //set to false to indicate that it is already compiled
				
				this.isCompiledToStaticFiled = true;
			} else {
				this.bytecodeGen(mv, argsMap, 2, refsMap, genClassName); //2 for args: double apply(Element e, Node n, double ...args);
			}
			mv.visitInsn(retType.getOpcode(Opcodes.IRETURN));
			
			mv.visitLocalVariable("this", "L" + genClassName + ";", null, startMatchesLabel, endMatchesLabel, 0);
			mv.visitLocalVariable("ap", param1.getDescriptor(),     null, startMatchesLabel, endMatchesLabel, 1);
			mv.visitLocalVariable("args", param2.getDescriptor(),   null, startMatchesLabel, endMatchesLabel, 2);
			mv.visitMaxs(-1, -1); // Auto generated
			cgen.endCode(mv, endMatchesLabel);

			cgen.endClass();

			byte[] bcode = cgen.dump();
			if (writeFile) {
				File dir = new File("gen_classes");
				if(!dir.exists()) {
					dir.mkdir();
				}

				FileOutputStream fos = new FileOutputStream("gen_classes/"+genClassName + ".class");
				fos.write(bcode);
				fos.close();
			}

			Class<?> c = mcl.defineClassForName(null, bcode);

			CompiledFunc func = (CompiledFunc) c.newInstance();

			// Set funcRefs field in CompiledFunc
			List<MathFunc> list = new ArrayList<MathFunc>();
			BytecodeUtils.postOrder(this, list);
			func.setFuncRefs(list.toArray(new MathFunc[0]));

			return func;
		} catch (Exception e) {
			e.printStackTrace();
			throw new RuntimeException(e.getMessage());
		}
	}

	// Indicate whether the result of the expression is needed to be assigned to a static field
	protected boolean compileToStaticField = false;
	// This flag is automatically set to false when 'compileToStaticField' is set to true and
	// it is set to true after the first time of compilation so the expressions that contain this
	// expression will refer the static field thereafter
	protected boolean isCompiledToStaticFiled = false;
	protected String genClassName;
	protected String staticFieldName;
	
	@Override
	public void compileToStaticField(boolean flag) {
		this.compileToStaticField = flag;
		this.isCompiledToStaticFiled = false;
	}

	//////////////Operator overloading support through Java-OO//////////////////

	public MathFunc valueOf(int v) {
		return new FC(v);
	}
	public MathFunc valueOf(long v) {
		return new FC(v);
	}
	public MathFunc valueOf(float v) {
		return new FC(v);
	}
	public MathFunc valueOf(double v) {
		return new FC(v);
	}

	public MathFunc add(MathFunc other) {
		return this.A(other);
	}
	public MathFunc add(int other) {
		return this.A(other);
	}
	public MathFunc addRev(int other) {
		return this.A(other);
	}
	public MathFunc add(long other) {
		return this.A(other);
	}
	public MathFunc addRev(long other) {
		return this.A(other);
	}	
	public MathFunc add(float other) {
		return this.A(other);
	}
	public MathFunc addRev(float other) {
		return this.A(other);
	}	
	public MathFunc add(double other) {
		return this.A(other);
	}
	public MathFunc addRev(double other) {
		return this.A(other);
	}
	
	public MathFunc subtract(MathFunc other) {
		return this.S(other);
	}
	public MathFunc subtract(int other) {
		return this.S(other);
	}
	public MathFunc subtractRev(int other) {
		return new FC(other).S(this);
	}
	public MathFunc subtract(long other) {
		return this.S(other);
	}
	public MathFunc subtractRev(long other) {
		return new FC(other).S(this);
	}	
	public MathFunc subtract(float other) {
		return this.S(other);
	}
	public MathFunc subtractRev(float other) {
		return new FC(other).S(this);
	}
	public MathFunc subtract(double other) {
		return this.S(other);
	}
	public MathFunc subtractRev(double other) {
		return new FC(other).S(this);
	}
	
	public MathFunc multiply(MathFunc other) {
		return this.M(other);
	}
	public MathFunc multiply(int other) {
		return this.M(other);
	}
	public MathFunc multiplyRev(int other) {
		return this.M(other);
	}
	public MathFunc multiply(long other) {
		return this.M(other);
	}
	public MathFunc multiplyRev(long other) {
		return this.M(other);
	}
	public MathFunc multiply(float other) {
		return this.M(other);
	}
	public MathFunc multiplyRev(float other) {
		return this.M(other);
	}
	public MathFunc multiply(double other) {
		return this.M(other);
	}
	public MathFunc multiplyRev(double other) {
		return this.M(other);
	}
	
	public MathFunc divide(MathFunc other) {
		return this.D(other);
	}	
	public MathFunc divide(int other) {
		return this.D(other);
	}
	public MathFunc divideRev(int other) {
		return new FC(other).D(this);
	}
	public MathFunc divide(long other) {
		return this.D(other);
	}
	public MathFunc divideRev(long other) {
		return new FC(other).D(this);
	}
	public MathFunc divide(float other) {
		return this.D(other);
	}
	public MathFunc divideRev(float other) {
		return new FC(other).D(this);
	}
	public MathFunc divide(double other) {
		return this.D(other);
	}
	public MathFunc divideRev(double other) {
		return new FC(other).D(this);
	}
	
	public MathFunc negate() {
		return new FSub(FMath.C0, this);
	};
	
	/////////////////////////////////////////////////////////////
	@Override
	public String getName() {
		return "[Undefined expressin. Please override 'String getName()']";
	}
	
	@Override
	public String getExpr() {
		return "[Undefined expressin. Please override 'String getExpr()']";
	}
	
	@Override
	public int getOpOrder() {
		return OP_ORDER0;
	}

	@Override
	public void setOpOrder(int order) {
		throw new UnsupportedOperationException();
	}
	
	@Override
	public String toString() {
		List<String> vars = this.getVarNames();
		if(vars.size() == 0) {
			if(getName() == null || getName().length() == 0)
				return "f() = " + getExpr();
			else
				return getName() + "() = " + getExpr();
		} else {
			StringBuilder sb = new  StringBuilder();
			sb.append("(");
			for(String arg : vars) {
				sb.append(arg).append(",");
			}
			sb.delete(sb.length()-1, sb.length());
			sb.append(")");
			if(getName() == null || getName().length() == 0)
				return "f" + sb.toString() + " = " + getExpr();
			else
				return getName() + sb.toString() + " = " + getExpr();
		}
	}
	
	@Override
	public MathFunc copy() {
		try {
			return (MathFunc)this.clone();
		} catch (CloneNotSupportedException e) {
			e.printStackTrace();
		}
		return null;
	}
	
	@Override
	public MathFunc setActiveVarByNames(List<String> varNames) {
		throw new UnsupportedOperationException();
	}
	
	@Override
	public List<String> getActiveVarNames() {
		throw new UnsupportedOperationException();
	}
	
	@Override
	public MathFunc setOuterVarActive() {
		throw new UnsupportedOperationException();
	}
	
	@Override
	public MathFunc setInnerVarActive() {
		throw new UnsupportedOperationException();
	}
	
	@Override
	public boolean isOuterVarActive() {
		throw new UnsupportedOperationException();
	}
	
	@Override
	public boolean isInnerVarActive() {
		throw new UnsupportedOperationException();
	}
}