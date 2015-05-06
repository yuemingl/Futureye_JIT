package edu.uta.futureye.function.basic;

import java.util.HashMap;
import java.util.Map;

import com.sun.org.apache.bcel.internal.generic.ConstantPoolGen;
import com.sun.org.apache.bcel.internal.generic.InstructionFactory;
import com.sun.org.apache.bcel.internal.generic.InstructionHandle;
import com.sun.org.apache.bcel.internal.generic.InstructionList;
import com.sun.org.apache.bcel.internal.generic.MethodGen;
import com.sun.org.apache.bcel.internal.generic.PUSH;

import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.function.AbstractMathFun;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.VariableArray;
import edu.uta.futureye.function.intf.MathFunc;

/**
 * Constant function: f = c
 * 
 */
public class FC extends AbstractMathFun{
	//Predefined constant
	public static FC C0 = new FC(0.0);
	public static FC C1 = new FC(1.0);
	public static FC Cm1 = new FC(-1.0);
	public static FC PI = new FC(Math.PI);
	public static FC E = new FC(Math.E);
	
	protected double val = 0.0;
	protected static Map<Double, FC> cs = new HashMap<Double, FC>();
	
	public FC() {
	}
	
	public FC(double val) {
		this.val = val;
	}
	
	/**
	 * 返回常值函数，并且保存在静态Map中，以便多次使用，节省内存
	 * 注意：不要使用该函数生成大量常数，否则内存在程序退出前不会释放
	 * @param v
	 * @return
	 */
	public static FC c(double v) {
		FC c = cs.get(v);
		if(c == null) {
			c = new FC(v);
			cs.put(v, c);
			return c;
		} else {
			return c;
		}
	}
	
	@Override
	public double apply() {
		return val;
	}
	
	@Override
	public double apply(Variable v) {
		return val;
	}

	@Override
	public double apply(Element e, Node n, double... args) {
		return val;
	}

	@Override
	public double apply(double... args) {
		return apply(null, null, args);
	}
	
	@Override
	public double apply(Variable v, Map<Object,Object> cache) {
		return val;
	}
	
	@Override
	public double[] applyAll(VariableArray v, Map<Object,Object> cache) {
		int len = v.length();
		double[] rlt = new double[len];
		for(int i=0;i<len;i++)
			rlt[i] = val;
		return rlt;
	}
	
	@Override
	public MathFunc _d(String varName) {
		return C0;
	}
	
	@Override
	public int getOpOrder() {
		return OP_ORDER0;
	}
	
	@Override
	public MathFunc copy() {
		return new FC(this.val);
	}
	
	@Override
	public String toString() {
		return String.valueOf(val);
	}
	
	@Override 
	public boolean isConstant() {
		return true;
	}
	
	@Override
	public InstructionHandle bytecodeGen(String clsName, MethodGen mg,
			ConstantPoolGen cp, InstructionFactory factory,
			InstructionList il, Map<String, Integer> argsMap, int argsStartPos, 
			Map<MathFunc, Integer> funcRefsMap) {
		return il.append(new PUSH(cp, val));
	}
}
