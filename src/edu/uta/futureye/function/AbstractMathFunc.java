package edu.uta.futureye.function;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import static com.sun.org.apache.bcel.internal.Constants.*;

import com.sun.org.apache.bcel.internal.generic.AALOAD;
import com.sun.org.apache.bcel.internal.generic.ALOAD;
import com.sun.org.apache.bcel.internal.generic.ArrayType;
import com.sun.org.apache.bcel.internal.generic.ClassGen;
import com.sun.org.apache.bcel.internal.generic.ConstantPoolGen;
import com.sun.org.apache.bcel.internal.generic.FieldGen;
import com.sun.org.apache.bcel.internal.generic.GETFIELD;
import com.sun.org.apache.bcel.internal.generic.InstructionFactory;
import com.sun.org.apache.bcel.internal.generic.InstructionHandle;
import com.sun.org.apache.bcel.internal.generic.InstructionList;
import com.sun.org.apache.bcel.internal.generic.MethodGen;
import com.sun.org.apache.bcel.internal.generic.PUSH;
import com.sun.org.apache.bcel.internal.generic.Type;
import com.sun.org.apache.bcel.internal.Constants;

import edu.uta.futureye.bytecode.CompiledFunc;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.util.BytecodeUtils;
import edu.uta.futureye.util.Constant;
import edu.uta.futureye.util.FuncClassLoader;

/**
 * Template to implement multiple variable MathFunc 
 *
 */
public abstract class AbstractMathFunc extends MathFuncBasic {
	protected List<String> varNames = new LinkedList<String>();
	protected int[] argIdx;
	protected String fName = null;
	
	public AbstractMathFunc() {
	}
	
	public AbstractMathFunc(List<String> varNames) {
		this.varNames = varNames;
	}
	
	public AbstractMathFunc(String varName, String ...aryVarNames) {
		varNames.add(varName);
		for(String s : aryVarNames)
			varNames.add(s);
	}
	
	@Override
	public List<String> getVarNames() {
		return varNames;
	}

	@Override
	public MathFunc setVarNames(List<String> varNames) {
		this.varNames = varNames;
		return this;
	}

	public MathFunc setArgIdx(int ...argIdx) {
		this.argIdx = argIdx;
		return this;
	}

	@Override
	public String getName() {
		return this.fName;
	}
	
	@Override
	public MathFunc setName(String name) {
		this.fName = name;
		return this;
	}
	
	@Override
	public boolean isConstant() {
		return false;
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
	public String getExpr() {
		String varList = getVarNames().toString();
		String displayVarList = "("+varList.substring(1, varList.length()-1)+")";
		
		Class<?> enclosingClass = getClass().getEnclosingClass();
		if (enclosingClass != null) {
		  return enclosingClass.getSimpleName() + displayVarList;
		} else {
		  return getClass().getSimpleName() + displayVarList;
		}
	}
	
	@Override
	public String toString() {
		if(getName() == null) {
			return getExpr();
		} else 
			return getName();
	}
}
