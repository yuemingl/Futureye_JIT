package edu.uta.futureye.util;

import static com.sun.org.apache.bcel.internal.Constants.*;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import com.sun.org.apache.bcel.internal.generic.ArrayType;
import com.sun.org.apache.bcel.internal.generic.ClassGen;
import com.sun.org.apache.bcel.internal.generic.ConstantPoolGen;
import com.sun.org.apache.bcel.internal.generic.FieldGen;
import com.sun.org.apache.bcel.internal.generic.InstructionConstants;
import com.sun.org.apache.bcel.internal.generic.InstructionFactory;
import com.sun.org.apache.bcel.internal.generic.InstructionList;
import com.sun.org.apache.bcel.internal.generic.MethodGen;
import com.sun.org.apache.bcel.internal.generic.PUTFIELD;
import com.sun.org.apache.bcel.internal.generic.Type;
import com.sun.org.apache.xpath.internal.operations.Variable;

import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.function.FBinaryOp;
import edu.uta.futureye.function.FCompose;
import edu.uta.futureye.function.intf.MathFunc;

public class BytecodeUtils {
	public static void postOrder(MathFunc func, List<MathFunc> list) {
		if(func instanceof FBinaryOp) {
			postOrder(((FBinaryOp) func).arg1, list);
			postOrder(((FBinaryOp) func).arg2, list);
		} else if(func instanceof FCompose) {
			FCompose fc = (FCompose)func;
			for(Entry<String, MathFunc> e : fc.fInners.entrySet()) {
				postOrder(e.getValue(), list);
			}
		}
		list.add(func);
	}
	
	public static Map<MathFunc, Integer> getFuncRefsMap(MathFunc func) {
		List<MathFunc> list = new ArrayList<MathFunc>();
		postOrder(func, list);
		Map<MathFunc, Integer> map = new HashMap<MathFunc, Integer>();
		for(int i=0; i<list.size(); i++) {
			map.put(list.get(i), i);
		}
		return map;
	}
	
	public static ClassGen genClass(MathFunc func, String[] varNames, String funcClsName, 
			boolean writeClassFile, boolean staticMethod) {
		
		String packageName = "edu.uta.futureye.bytecode";
		String clsName = funcClsName;
		String fullClsName = packageName+"."+clsName;
		ClassGen cg = new ClassGen(fullClsName, "edu.uta.futureye.bytecode.CompiledFunc",
				"<generated>", ACC_PUBLIC | ACC_SUPER, null);
		ConstantPoolGen cp = cg.getConstantPool(); // cg creates constant pool
		InstructionList il = new InstructionList();
		InstructionFactory factory = new InstructionFactory(cg);
		

		short acc_flags = ACC_PUBLIC;
		if(staticMethod)
			acc_flags |= ACC_STATIC;
		MethodGen mg = new MethodGen(acc_flags, // access flags
				Type.DOUBLE, // return type
				new Type[] { // argument types
					Type.getType(Element.class),
					Type.getType(Node.class),
					new ArrayType(Type.DOUBLE, 1) 
				}, 
				new String[] { "e", "n", "args" }, // arg names
				"apply", fullClsName, // method, class
				il, cp);
		
		System.out.println("JIT Compiled: "+func);

		
		HashMap<String, Integer> argsMap = new HashMap<String, Integer>();
		if(varNames == null) {
			List<String> args = func.getVarNames();
			for(int i=0; i<args.size(); i++) {
				argsMap.put(args[i], i);
			}
		} else {
			for(int i=0; i<varNames.length; i++) {
				argsMap.put(varNames[i], i);
			}
			func.setArgIdx(argsMap);
		}


		Map<MathFunc, Integer> refsMap = getFuncRefsMap(func);
		
		if(staticMethod)
			func.bytecodeGen(clsName, mg, cp, factory, il, argsMap, 2, refsMap);
		else
			func.bytecodeGen(clsName, mg, cp, factory, il, argsMap, 3, refsMap);
		il.append(InstructionConstants.DRETURN);
		
//	Test
//		FieldGen fg = new FieldGen(ACC_PUBLIC, new ArrayType(Type.getType(MathFunc.class), 1), "funcRefs", cp);
//		il.append(InstructionConstants.ALOAD_0);
//		il.append(InstructionConstants.ACONST_NULL);
//		il.append(new PUTFIELD(cp.addFieldref(clsName, "funcRefs", fg.getSignature())));
//		cg.addField(fg.getField());

//		InstructionList il = new InstructionList();
//		il.append(new ALOAD(0));
//		il.append(new ALOAD(0));
//		il.append(new GETFIELD(cp.addFieldref(name, "_counter", "I")));
//		il.append(new ICONST(1));
//		il.append(new IADD());
//		il.append(new PUTFIELD(cp.addFieldref(name, "_counter", "I")));
//		mg.getInstructionList().insert(il);
		
		
		mg.setMaxStack();
		cg.addMethod(mg.getMethod());
		il.dispose(); // Allow instruction handles to be reused
		
		cg.addEmptyConstructor(ACC_PUBLIC);
		if(writeClassFile) {
			try {
				cg.getJavaClass().dump("bin/edu/uta/futureye/bytecode/"+clsName+".class");
			} catch (java.io.IOException e) {
				System.err.println(e);
			}
		}
		
		
		return cg;
	}
}
