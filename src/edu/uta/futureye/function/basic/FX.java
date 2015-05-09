package edu.uta.futureye.function.basic;

import static com.sun.org.apache.bcel.internal.generic.InstructionConstants.DALOAD;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
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
import edu.uta.futureye.function.MathFuncBasic;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.VariableArray;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.util.Constant;

/**
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
	int argIdx;
	
	/**
	 * Identity function: f(varName) = varName
	 */
	public FX(String varName) {
		this.varName = varName;
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
	public InstructionHandle bytecodeGen(String clsName, MethodGen mg,
			ConstantPoolGen cp, InstructionFactory factory,
			InstructionList il, Map<String, Integer> argsMap, int argsStartPos, 
			Map<MathFunc, Integer> funcRefsMap) {
		il.append(new ALOAD(argsStartPos));
		il.append(new PUSH(cp, argsMap.get(varName)));
		return il.append(DALOAD);
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
}
