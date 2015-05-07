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

public abstract class AbstractMathFun implements MathFunc {
	protected List<String> varNames = new LinkedList<String>();
	protected String fName = null;
	
	public AbstractMathFun() {
	}
	
	public AbstractMathFun(List<String> varNames) {
		this.varNames = varNames;
	}
	
	public AbstractMathFun(String varName, String ...aryVarNames) {
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
	
	/**
	 * Implement this function yourself
	 */
	@Override
	public abstract double apply(Variable v);
	
	public abstract double apply(double ...args);
	
	public double apply(Element e, Node n, double... args) {
		return apply(args);
	}
	
	@Override
	public double apply(Variable v, Map<Object,Object> cache) {
		//Ignore cache
		return apply(v);
	}
	
	@Override
	public double[] applyAll(VariableArray v, Map<Object,Object> cache) {
		throw new UnsupportedOperationException();
	}

	@Override
	public double apply() {
		throw new UnsupportedOperationException();
	}
	
	@Override
	public double apply(double x) {
		throw new UnsupportedOperationException();
	}
	
	@Override
	public double apply(double x, double y) {
		throw new UnsupportedOperationException();
	}
	
	@Override
	public double apply(double x, double y, double z) {
		throw new UnsupportedOperationException();
	}
	
	/**
	 * Implement this function yourself if necessary
	 */
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
		else
			return new FCompose(this, fInners);
	}
	
	////////////////////////Operations////////////////////////////////////
	
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
			return FC.C0;
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
			return FC.C0;
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
	
	@Override
	public MathFunc copy() {
		throw new UnsupportedOperationException();
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
	public int getOpOrder() {
		return OP_ORDER3;
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
	
	@Override
	public boolean isConstant() {
		return false;
	}
	
	///////////////////////////Compilation///////////////////////////////
	
	/**
	 * 
	 */
	@Override
	public InstructionHandle bytecodeGen(String clsName, MethodGen mg, 
			ConstantPoolGen cp, InstructionFactory factory, 
			InstructionList il, Map<String, Integer> argsMap, 
			int argsStartPos, Map<MathFunc, Integer> funcRefsMap) {

		FieldGen fg = new FieldGen(ACC_PUBLIC, new ArrayType(Type.getType(MathFunc.class), 1), "funcRefs", cp);
		//System.out.println(fg.getSignature());
		int idxFuncRefs = cp.addFieldref("edu.uta.futureye.bytecode.CompiledFunc", "funcRefs", fg.getSignature());
		
		il.append(new ALOAD(0));
		il.append(new GETFIELD(idxFuncRefs));
		il.append(new PUSH(cp, funcRefsMap.get(this)));
		il.append(new AALOAD());
		il.append(new ALOAD(1));
		il.append(new ALOAD(2));
		il.append(new ALOAD(3));
		return  il.append(factory.createInvoke("edu.uta.futureye.function.intf.MathFunc", "apply",
				Type.DOUBLE, new Type[] {
					Type.getType(Element.class),
					Type.getType(Node.class),
					new ArrayType(Type.DOUBLE, 1)
				}, 
				Constants.INVOKEINTERFACE));
	}
	
	/**
	 * 
	 */
	@Override
	public CompiledFunc compile() {
		String clsName = this.fName + java.util.UUID.randomUUID().toString().replaceAll("-", "");
		
		FuncClassLoader<CompiledFunc> fcl = new FuncClassLoader<CompiledFunc>();
		ClassGen genClass = BytecodeUtils.genClass(this, clsName, true, false);
		CompiledFunc func = fcl.newInstance(genClass);
		
		List<MathFunc> list = new ArrayList<MathFunc>();
		BytecodeUtils.postOrder(this, list);
		func.setFuncRefs(list.toArray(new MathFunc[0]));
		
		return func;
	}
	
	//////////////Operator overloading support through Java-OO//////////////////
	/**
	 * Operator overloading support:
	 * MathFunc a = 5;
	 * 
	 * @param v
	 * @return
	 */
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
	/**
	 * Operator overload support:
	 * a+b
	 * @param other
	 * @return
	 */
	public MathFunc add(MathFunc other) {
		return new FAdd(this, other);
	}
	public MathFunc add(int other) {
		return new FAdd(this, new FC(other));
	}
	public MathFunc addRev(int other) {
		return new FAdd(new FC(other), this);
	}
	public MathFunc add(long other) {
		return new FAdd(this, new FC(other));
	}
	public MathFunc addRev(long other) {
		return new FAdd(new FC(other), this);
	}	
	public MathFunc add(float other) {
		return new FAdd(this, new FC(other));
	}
	public MathFunc addRev(float other) {
		return new FAdd(new FC(other), this);
	}	
	public MathFunc add(double other) {
		return new FAdd(this, new FC(other));
	}
	public MathFunc addRev(double other) {
		return new FAdd(new FC(other), this);
	}
	
	/**
	 * Operator overload support:
	 * a-b
	 * @param other
	 * @return
	 */
	public MathFunc subtract(MathFunc other) {
		return new FSub(this, other);
	}
	public MathFunc subtract(int other) {
		return new FSub(this, new FC(other));
	}
	public MathFunc subtractRev(int other) {
		return new FSub(new FC(other), this);
	}
	public MathFunc subtract(long other) {
		return new FSub(this, new FC(other));
	}
	public MathFunc subtractRev(long other) {
		return new FSub(new FC(other), this);
	}	
	public MathFunc subtract(float other) {
		return new FSub(this, new FC(other));
	}
	public MathFunc subtractRev(float other) {
		return new FSub(new FC(other), this);
	}
	public MathFunc subtract(double other) {
		return new FSub(this, new FC(other));
	}
	public MathFunc subtractRev(double other) {
		return new FSub(new FC(other), this);
	}
	
	/**
	 * Operator overload support:
	 * a*b
	 * @param other
	 * @return
	 */
	public MathFunc multiply(MathFunc other) {
		return new FMul(this, other);
	}
	public MathFunc multiply(int other) {
		return new FMul(this, new FC(other));
	}
	public MathFunc multiplyRev(int other) {
		return new FMul(new FC(other), this);
	}
	public MathFunc multiply(long other) {
		return new FMul(this, new FC(other));
	}
	public MathFunc multiplyRev(long other) {
		return new FMul(new FC(other), this);
	}
	public MathFunc multiply(float other) {
		return new FMul(this, new FC(other));
	}
	public MathFunc multiplyRev(float other) {
		return new FMul(new FC(other), this);
	}
	public MathFunc multiply(double other) {
		return new FMul(this, new FC(other));
	}
	public MathFunc multiplyRev(double other) {
		return new FMul(new FC(other), this);
	}
	
	/**
	 * Operator overload support:
	 * a/b
	 * @param other
	 * @return
	 */
	public MathFunc divide(MathFunc other) {
		return new FDiv(this, other);
	}	
	public MathFunc divide(int other) {
		return new FDiv(this, new FC(other));
	}
	public MathFunc divideRev(int other) {
		return new FDiv(new FC(other), this);
	}
	public MathFunc divide(long other) {
		return new FDiv(this, new FC(other));
	}
	public MathFunc divideRev(long other) {
		return new FDiv(new FC(other), this);
	}
	public MathFunc divide(float other) {
		return new FDiv(this, new FC(other));
	}
	public MathFunc divideRev(float other) {
		return new FDiv(new FC(other), this);
	}
	public MathFunc divide(double other) {
		return new FDiv(this, new FC(other));
	}
	public MathFunc divideRev(double other) {
		return new FDiv(new FC(other), this);
	}
	
	/**
	 * Operator overload support:
	 * -a
	 * 
	 */
	public MathFunc negate() {
		return new FSub(FC.C0, this);
	};

}
