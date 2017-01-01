package edu.uta.futureye.test;

import java.util.Map;

import org.apache.bcel.generic.ALOAD;
import org.apache.bcel.generic.ConstantPoolGen;
import org.apache.bcel.generic.DALOAD;
import org.apache.bcel.generic.DMUL;
import org.apache.bcel.generic.InstructionFactory;
import org.apache.bcel.generic.InstructionHandle;
import org.apache.bcel.generic.InstructionList;
import org.apache.bcel.generic.MethodGen;
import org.apache.bcel.generic.PUSH;

import edu.uta.futureye.bytecode.CompiledFunc;
import edu.uta.futureye.function.SingleVarFunc;
import edu.uta.futureye.function.intf.MathFunc;

public class TestSetArgIdxAndCodeGen {

	/**
	 * f(x)=x*x
	 * @author yueming.liu
	 *
	 */
	public static class MyFunc1 extends SingleVarFunc {
		public MyFunc1(String funcName, String varName) {
			super(funcName, varName);
		}
		@Override
		public double apply(double... args) {
			System.out.println("MyFunc1.apply(): this is called");
			return args[this.argIdx] * args[this.argIdx];
		}
		
		@Override
		public String getExpr() {
			return "x*x";
		}
		
	}
	
	/**
	 * f(x)=x*x
	 * compiled version
	 * @author yueming.liu
	 *
	 */
	public static class MyFunc2 extends SingleVarFunc {
		public MyFunc2(String funcName, String varName) {
			super(funcName, varName);
		}
		@Override
		public double apply(double... args) {
			System.out.println("MyFunc2.apply(): this function should not be called");
			return args[this.argIdx] * args[this.argIdx];
		}
		
		/**
		 * generate byte code for f(x)=x*x
		 * the function 'apply()' will not be call by override this function
		 */
		@Override
		public InstructionHandle bytecodeGen(String clsName, MethodGen mg, 
				ConstantPoolGen cp, InstructionFactory factory, 
				InstructionList il, Map<String, Integer> argsMap, 
				int argsStartPos, Map<MathFunc, Integer> funcRefsMap) {
			il.append(new ALOAD(argsStartPos));
			il.append(new PUSH(cp, argsMap.get(this.getVarName())));
			il.append(new DALOAD());
			il.append(new ALOAD(argsStartPos));
			il.append(new PUSH(cp, argsMap.get(this.getVarName())));
			il.append(new DALOAD());
			return il.append(new DMUL());
		}
		
		@Override
		public String getExpr() {
			return "x*x";
		}
	}

	public static void main(String[] args) {
		MyFunc1 func1 = new MyFunc1("func1","x");
		System.out.println(func1);
		System.out.println(func1.apply(new double[]{3}));
		CompiledFunc cfunc1 = func1.compile();
		System.out.println(cfunc1.apply(new double[]{3}));
		
		MyFunc2 func2 = new MyFunc2("func2","x");
		System.out.println(func2);
		CompiledFunc cfunc2 = func2.compile();
		System.out.println(cfunc2.apply(new double[]{3}));
		
	}

}
