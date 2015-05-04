package edu.uta.futureye.function;

import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import com.sun.org.apache.bcel.internal.generic.ConstantPoolGen;
import com.sun.org.apache.bcel.internal.generic.InstructionFactory;
import com.sun.org.apache.bcel.internal.generic.InstructionList;
import com.sun.org.apache.bcel.internal.generic.MethodGen;

import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.intf.MathFun;
import edu.uta.futureye.util.Constant;

public abstract class AbstractMathFun implements MathFun {
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
	public MathFun setVarNames(List<String> varNames) {
		this.varNames = varNames;
		return this;
	}
	
	/**
	 * Implement this function yourself
	 */
	@Override
	public abstract double apply(Variable v);
	
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
	public MathFun _d(String varName) {
		throw new UnsupportedOperationException();
	}

	@Override
	public MathFun compose(final Map<String,MathFun> fInners) {
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
	public MathFun A(MathFun g) {
		final MathFun f1 = this;
		final MathFun f2 = g;
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
	public MathFun A(double g) {
		return A(FC.c(g));
	}
	
	@Override
	public MathFun S(MathFun g) {
		final MathFun f1 = this;
		final MathFun f2 = g;
		if(f1.isConstant() && f2.isConstant()) {
			return new FC(f1.apply() - f2.apply());
		} else if(f2.isConstant() && Math.abs(f2.apply()) < Constant.eps) {
			return f1;
		} else {
			return new FSub(f1, f2);
		}
	}
	@Override
	public MathFun S(double g) {
		return S(FC.c(g));
	}	
	
	@Override
	public MathFun M(MathFun f) {
		final MathFun f1 = this;
		final MathFun f2 = f;
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
	public MathFun M(double g) {
		return M(FC.c(g));
	}	
	
	@Override
	public MathFun D(MathFun f) {
		final MathFun f1 = this;
		final MathFun f2 = f;
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
	public MathFun D(double g) {
		return D(FC.c(g));
	}
	
	@Override
	public MathFun copy() {
		throw new UnsupportedOperationException();
	}
	
	@Override
	public String getFunName() {
		return this.fName;
	}
	
	@Override
	public MathFun setFunName(String name) {
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
	public String getExpression() {
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
		if(getFunName() == null) {
			return getExpression();
		} else 
			return getFunName();
	}
	
	@Override
	public boolean isConstant() {
		return false;
	}

}
