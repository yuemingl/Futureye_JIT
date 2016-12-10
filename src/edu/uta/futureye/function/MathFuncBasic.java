package edu.uta.futureye.function;

import static org.apache.bcel.Constants.ACC_PUBLIC;

import java.util.ArrayList;
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

import edu.uta.futureye.bytecode.CompiledFunc;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.basic.FComposite;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.function.operator.FAdd;
import edu.uta.futureye.function.operator.FDiv;
import edu.uta.futureye.function.operator.FMul;
import edu.uta.futureye.function.operator.FSub;
import edu.uta.futureye.util.BytecodeUtils;
import edu.uta.futureye.util.Constant;
import edu.uta.futureye.util.FuncClassLoader;

public abstract class MathFuncBasic implements MathFunc, Cloneable {
	@Override
	public abstract double apply(double ...args);
	
	@Override
	public double apply(Element e, Node n, double... args) {
		return apply(args);
	}

	@Deprecated
	@Override
	public double apply(Variable v) {
		Node n = new Node(v.getIndex());
		return apply(v.getElement(), n, v.getVarValues());
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
		else
			return new FComposite(this, fInners);
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
	 * Call the 'apply' method by default if a sub-class does not override this method to generate the bytecode
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
	
	@Override
	public CompiledFunc compile() {
		String clsName="";
		if(getName() == null || getName().length() == 0)
			clsName = this.getClass().getSimpleName();
		clsName = clsName + java.util.UUID.randomUUID().toString().replaceAll("-", "");
		
		FuncClassLoader<CompiledFunc> fcl = new FuncClassLoader<CompiledFunc>();
		ClassGen genClass = BytecodeUtils.genClass(this, null, clsName, true, false);
		CompiledFunc func = fcl.newInstance(genClass);
		
		List<MathFunc> list = new ArrayList<MathFunc>();
		BytecodeUtils.postOrder(this, list);
		func.setFuncRefs(list.toArray(new MathFunc[0]));
		
		return func;
	}
	
	@Override
	public CompiledFunc compile(String[] varNames) {
		String clsName = getName();
		if(clsName == null || clsName.length() == 0)
			clsName = this.getClass().getSimpleName();
		clsName = clsName + java.util.UUID.randomUUID().toString().replaceAll("-", "");
		
		FuncClassLoader<CompiledFunc> fcl = new FuncClassLoader<CompiledFunc>();
		ClassGen genClass = BytecodeUtils.genClass(this, varNames, clsName, true, false);
		CompiledFunc func = fcl.newInstance(genClass);
		
		List<MathFunc> list = new ArrayList<MathFunc>();
		BytecodeUtils.postOrder(this, list);
		func.setFuncRefs(list.toArray(new MathFunc[0]));
		
		return func;
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
		return "";
	}
	
	@Override
	public String getExpr() {
		return "";
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
			if(getName().length() == 0)
				return getExpr();
			else
				return getName() + " = " + getExpr();
		} else {
			StringBuilder sb = new  StringBuilder();
			sb.append("(");
			for(String arg : vars) {
				sb.append(arg).append(",");
			}
			sb.delete(sb.length()-1, sb.length());
			sb.append(")");
			if(getName().length() == 0)
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
	public MathFunc setActiveVarNames(List<String> varNames) {
		//Do nothing
		return this;
	}

}