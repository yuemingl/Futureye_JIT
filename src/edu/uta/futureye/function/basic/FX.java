package edu.uta.futureye.function.basic;

import static com.sun.org.apache.bcel.internal.generic.InstructionConstants.DALOAD;

import java.util.Map;

import com.sun.org.apache.bcel.internal.generic.ALOAD;
import com.sun.org.apache.bcel.internal.generic.ConstantPoolGen;
import com.sun.org.apache.bcel.internal.generic.InstructionFactory;
import com.sun.org.apache.bcel.internal.generic.InstructionHandle;
import com.sun.org.apache.bcel.internal.generic.InstructionList;
import com.sun.org.apache.bcel.internal.generic.MethodGen;
import com.sun.org.apache.bcel.internal.generic.PUSH;

import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.function.AbstractSimpleMathFunc;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.VariableArray;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.util.Constant;

/**
 * f(x) = x
 * 
 */
public class FX extends AbstractSimpleMathFunc {
	/**
	 * Predefined instances of FX
	 */
	public final static FX fx = new FX(Constant.x); 
	public final static FX fy = new FX(Constant.y); 
	public final static FX fz = new FX(Constant.z); 
	
	public final static FX fr = new FX(Constant.r); 
	public final static FX fs = new FX(Constant.s); 
	public final static FX ft = new FX(Constant.t); 
	
	/**
	 * Use this to construct a function: f(varName) = varName
	 */
	public FX(String varName) {
		super(varName, varName);
	}

	@Override
	public double apply(Variable v) {
		return v.get(varName);
	}

	@Override
	public double apply(Element e, Node n, double... args) {
		return args[argIdx];
	}

	@Override
	public double apply(double... args) {
		return apply(null, null, args);
	}
	
	@Override
	public double apply(Variable v, Map<Object,Object> cache) {
		return v.get(varName);
	}
	
	@Override
	public double[] applyAll(VariableArray v, Map<Object,Object> cache) {
		return v.get(varName);
	}

	@Override
	public MathFunc diff(String varName) {
		if(this.varName.equals(varName))
			return FC.C1;
		else
			return FC.C0;
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
	public MathFunc copy() {
		FX ret = new FX(this.varName);
		ret.fName = this.fName;
		ret.argIdx = this.argIdx;
		return ret;
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

}
