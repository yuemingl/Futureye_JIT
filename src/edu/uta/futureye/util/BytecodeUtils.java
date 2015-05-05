package edu.uta.futureye.util;

import static com.sun.org.apache.bcel.internal.Constants.ACC_PUBLIC;
import static com.sun.org.apache.bcel.internal.Constants.ACC_STATIC;
import static com.sun.org.apache.bcel.internal.Constants.ACC_SUPER;

import java.util.HashMap;
import java.util.List;

import com.sun.org.apache.bcel.internal.generic.ArrayType;
import com.sun.org.apache.bcel.internal.generic.ClassGen;
import com.sun.org.apache.bcel.internal.generic.ConstantPoolGen;
import com.sun.org.apache.bcel.internal.generic.InstructionConstants;
import com.sun.org.apache.bcel.internal.generic.InstructionFactory;
import com.sun.org.apache.bcel.internal.generic.InstructionList;
import com.sun.org.apache.bcel.internal.generic.MethodGen;
import com.sun.org.apache.bcel.internal.generic.Type;

import edu.uta.futureye.function.intf.MathFunc;

public class BytecodeUtils {
	public static ClassGen genClassBytecodeFunc(MathFunc func, String funcClsName, boolean writeClassFile, boolean staticMethod) {
		String packageName = "edu.uta.futureye.bytecode";
		String clsName = funcClsName;
		String fullClsName = packageName+"."+clsName;
		ClassGen cg = new ClassGen(fullClsName, "java.lang.Object",
				"<generated>", ACC_PUBLIC | ACC_SUPER, new String[]{"edu.uta.futureye.bytecode.BytecodeFunc"});
		ConstantPoolGen cp = cg.getConstantPool(); // cg creates constant pool
		InstructionList il = new InstructionList();
		InstructionFactory factory = new InstructionFactory(cg);
		
		short acc_flags = ACC_PUBLIC;
		if(staticMethod)
			acc_flags |= ACC_STATIC;
		MethodGen mg = new MethodGen(acc_flags, // access flags
				Type.DOUBLE, // return type
				new Type[] { // argument types
					new ArrayType(Type.DOUBLE, 1) 
				}, 
				new String[] { "args" }, // arg names
				"apply", fullClsName, // method, class
				il, cp);
		
		System.out.println("JIT Compiled: "+func);
		
		HashMap<String, Integer> argsMap = new HashMap<String, Integer>();
		List<String> args = func.getVarNames();
		for(int i=0; i<args.size(); i++) {
			argsMap.put(args[i], i);
		}
		if(staticMethod)
			func.bytecodeGen(mg, cp, factory, il, argsMap, 0);
		else
			func.bytecodeGen(mg, cp, factory, il, argsMap, 1);
		il.append(InstructionConstants.DRETURN);
		
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
