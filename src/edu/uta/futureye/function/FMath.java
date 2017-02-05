package edu.uta.futureye.function;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.apache.bcel.generic.ConstantPoolGen;
import org.apache.bcel.generic.InstructionFactory;
import org.apache.bcel.generic.InstructionHandle;
import org.apache.bcel.generic.InstructionList;
import org.apache.bcel.generic.MethodGen;

import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.basic.FX;
import edu.uta.futureye.function.basic.SpaceVectorFunction;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.function.intf.VectorMathFunc;
import edu.uta.futureye.function.operator.FAbs;
import edu.uta.futureye.function.operator.FAcos;
import edu.uta.futureye.function.operator.FAsin;
import edu.uta.futureye.function.operator.FCos;
import edu.uta.futureye.function.operator.FCosh;
import edu.uta.futureye.function.operator.FExp;
import edu.uta.futureye.function.operator.FLog;
import edu.uta.futureye.function.operator.FLog10;
import edu.uta.futureye.function.operator.FMax;
import edu.uta.futureye.function.operator.FMin;
import edu.uta.futureye.function.operator.FPow;
import edu.uta.futureye.function.operator.FSignum;
import edu.uta.futureye.function.operator.FSin;
import edu.uta.futureye.function.operator.FSinh;
import edu.uta.futureye.function.operator.FSqrt;
import edu.uta.futureye.function.operator.FTan;
import edu.uta.futureye.function.operator.FTanh;
import edu.uta.futureye.util.BytecodeUtils;
import edu.uta.futureye.util.Constant;
import edu.uta.futureye.util.FutureyeException;
import edu.uta.futureye.util.Utils;
import edu.uta.futureye.util.container.ObjList;

public class FMath {
	//--- Predefined static objects ------------------------
	/**
	 * Use "import static edu.uta.futureye.function.operator.FMath.*" to simplify
	 * the usage of these predefined static objects
	 */
	public final static FC C0 = new FC(0.0);
	public final static FC C1 = new FC(1.0);
	public final static FC Cm1 = new FC(-1.0);
	public final static FC PI = new FC(Math.PI);
	public final static FC E = new FC(Math.E);
	
	public final static FX x = new FX(Constant.x); 
	public final static FX y = new FX(Constant.y); 
	public final static FX z = new FX(Constant.z); 
	
	public final static FX r = new FX(Constant.r); 
	public final static FX s = new FX(Constant.s); 
	public final static FX t = new FX(Constant.t); 
	
	public static MathFunc C(double v) {
		return FC.c(v);
	}
	
	
	//--- Basic operations ------------------------
	public static MathFunc sin(MathFunc f) {
		return new FSin(f);
	}
	
	public static MathFunc cos(MathFunc f) {
		return new FCos(f);
	}
	
	public static MathFunc tan(MathFunc f) {
		return new FTan(f);
	}
	
	public static MathFunc sqrt(final MathFunc f) {
		return new FSqrt(f);
	}
	
	public static MathFunc pow(MathFunc base, double exp) {
		return new FPow(base, FC.c(exp));
	}
	
	public static MathFunc pow(MathFunc base, MathFunc exp) {
		return new FPow(base, exp);
	}
	
	public static MathFunc abs(MathFunc f) {
		return new FAbs(f);
	}
	
	public static MathFunc signum(MathFunc f) {
		return new FSignum(f);
	}
	
	public static MathFunc sinh(MathFunc f) {
		return new FSinh(f);
	}

	public static MathFunc cosh(MathFunc f) {
		return new FCosh(f);
	}

	public static MathFunc tanh(MathFunc f) {
		return new FTanh(f);
	}

	public static MathFunc asin(MathFunc f) {
		return new FAsin(f);
	}

	public static MathFunc acos(MathFunc f) {
		return new FAcos(f);
	}

	public static MathFunc log(MathFunc f) {
		return new FLog(f);
	}

	public static MathFunc log10(MathFunc f) {
		return new FLog10(f);
	}

	public static MathFunc max(MathFunc f,MathFunc g) {
		return new FMax(f,g);
	}

	public static MathFunc min(MathFunc f,MathFunc g) {
		return new FMin(f,g);
	}

	public static MathFunc exp(MathFunc f) {
		return new FExp(f);
	}
	
	/**
	 * \sum{f_i} = f_1 + f_2 + ... + f_N
	 * 
	 * @param fi
	 * @return
	 */
	public static MathFunc sum(MathFunc ...fi) {
		if(fi==null || fi.length==0) {
			throw new FutureyeException("Check parameter f="+fi);
		}
		MathFunc rlt = fi[0];
		for(int i=1;i<fi.length;i++) {
			rlt = rlt.A(fi[i]);
		}
		return rlt;
	}
	
	/**
	 * c1*f1 + c2*f2
	 * 
	 * @param c1
	 * @param f1
	 * @param c2
	 * @param f2
	 * @return
	 */
	public static MathFunc linearCombination(double c1, MathFunc f1,
			double c2,MathFunc f2) {
		return new FC(c1).M(f1).A(new FC(c2).M(f2));
	}
	
	static class FLinearCombination extends MultiVarFunc{
		double []ci;
		MathFunc []fi;
		public FLinearCombination(double []ci, MathFunc []fi) {
			int len = ci.length;
			this.ci = new double[len];
			this.fi = new MathFunc[len];
			List<String> list = new ArrayList<String>();
			for(int i=0;i<fi.length;i++) {
				this.ci[i] = ci[i];
				this.fi[i] = fi[i];
				list = Utils.mergeList(list, fi[i].getVarNames());
			}
			this.setVarNames(list);
			this.setArgIdx(Utils.getIndexMap(list));
		}
		
		@Override
		public double apply(Variable v) {
			double rlt = 0.0;
			for(int i=0;i<fi.length;i++) {
				rlt += ci[i]*fi[i].apply(v);
			}
			return rlt;
		}
		
		@Override
		public double apply(Variable v, Map<Object,Object> cache) {
			double rlt = 0.0;
			for(int i=0;i<fi.length;i++) {
				rlt += ci[i]*fi[i].apply(v,cache);
			}
			return rlt;
		}
		
		@Override
		public double[] applyAll(VariableArray v, Map<Object,Object> cache) {
			int len = v.length();
			double[] rlt = fi[0].applyAll(v,cache);
			for(int j=0;j<len;j++) 
				rlt[j] *= ci[0];
			for(int i=1;i<fi.length;i++) {
				double[] vs = fi[i].applyAll(v,cache);
				for(int j=0;j<len;j++) {
					rlt[j] += ci[i]*vs[j];
				}
			}
			return rlt;
		}
		
		@Override
		public MathFunc diff(String varName) {
			MathFunc[] fdi = new MathFunc[fi.length];
			for(int i=0;i<fi.length;i++) {
				//fdi[i] = fi[i].diff(varName).setVarNames(fi[i].getVarNames());
				fdi[i] = fi[i].diff(varName);
			}
			return new FLinearCombination(ci,fdi);
		}
		
		@Override
		public int getOpOrder() {
			return OP_ORDER3;
		}
		
		@Override
		public String getExpr() {
			StringBuilder sb = new StringBuilder();
			sb.append(ci[0]);
			sb.append("*");
			sb.append(fi[0].getExpr());
			for(int i=1;i<fi.length;i++) {
				sb.append(" + ");
				sb.append(ci[i]);
				sb.append("*");
				if(OP_ORDER2 < fi[i].getOpOrder())
					sb.append("(").append(fi[i].getExpr()).append(")");
				else
					sb.append(fi[i].getExpr());
			}
			return sb.toString();
		}
		
		@Override
		public String toString() {
			return getExpr();
		}

		@Override
		public double apply(double... args) {
			double rlt = 0.0;
			for(int i=0;i<fi.length;i++) {
				rlt += ci[i]*fi[i].apply(args);
			}
			return rlt;
		}
		
		@Override
		public InstructionHandle bytecodeGen(String clsName, MethodGen mg,
				ConstantPoolGen cp, InstructionFactory factory,
				InstructionList il, Map<String, Integer> argsMap, int argsStartPos, 
				Map<MathFunc, Integer> funcRefsMap) {
			MathFunc ret = C0;
			for(int i=0;i<fi.length;i++) {
				ret = ret + ci[i]*fi[i];
			}
			Map<MathFunc, Integer> refsMap2 = BytecodeUtils.getFuncRefsMap(ret);
			funcRefsMap.putAll(refsMap2);
			ret.setArgIdx(this.getArgIdxMap());
			return ret.bytecodeGen(clsName, mg, cp, factory, il, argsMap, argsStartPos, funcRefsMap);
		}
		
		@Override
		public MathFunc setArgIdx(Map<String, Integer> argsMap) {
			super.setArgIdx(argsMap);
			return this;
		}				
	}
	
	/**
	 * \sum{c_i*f_i} = c_1*f_1 + c_2*f_2 + ... + c_N*f_N
	 * 
	 * @param ci
	 * @param fi
	 * @return
	 */
	public static MathFunc linearCombination(double []ci, MathFunc []fi) {
		if(ci==null || ci.length==0 || fi==null || fi.length==0) {
			throw new FutureyeException("Check parameters ci="+ci+", fi="+fi);
		} else if(ci.length != fi.length) {
			throw new FutureyeException("(ci.length="+ci.length+") != (fi.lenght="+fi.length+")");
		}
		
//		Function rlt = new FC(ci[0]).M(fi[0]);
//		for(int i=1;i<fi.length;i++) {
//			rlt = rlt.A(new FC(ci[i]).M(fi[i]));
//		}
//		return rlt;

		return new FLinearCombination(ci,fi);
	}	
	
	
	public static MathFunc dot(VectorMathFunc f, VectorMathFunc g) {
		return f.dot(g);
	}
	
	
	/**
	 * Compute gradient of <code>fun</code>
	 * 
	 * @param fun
	 * @return
	 */
	public static VectorMathFunc grad(MathFunc fun) {
		List<String> names = fun.getVarNames();
		VectorMathFunc rlt = new SpaceVectorFunction(names.size());
		for(int i=0;i<names.size();i++)
			rlt.set(i+1, fun.diff(names.get(i)));
		return rlt;
	}
	
	/**
	 * Compute gradient of <code>fun</code> with respect to variables <code>varNames</code>
	 * in case of composition of functions.
	 * 
	 * @param vars
	 * @return
	 */
	public static VectorMathFunc grad(MathFunc fun, String ...varNames) {
		VectorMathFunc rlt = new SpaceVectorFunction(varNames.length);
		for(int i=0;i<varNames.length;i++)
			rlt.set(i+1, fun.diff(varNames[i]));
		return rlt;
	}
	
	/**
	 * Compute gradient of <code>fun</code> with respect to variables <code>varNames</code>
	 * in case of composition of functions.
	 * 
	 * @param fun
	 * @param varNames
	 * @return
	 */
	public static VectorMathFunc grad(MathFunc fun,ObjList<String> varNames) {
		VectorMathFunc rlt = new SpaceVectorFunction(varNames.size());
		for(int i=1;i<=varNames.size();i++)
			rlt.set(i, fun.diff(varNames.at(i)));
		return rlt;
	}
	
	/**
	 * Compute divergence of <code>vFun</code>
	 * 
	 * @param fun
	 * @return
	 */
	public static MathFunc div(VectorMathFunc vFun) {
		int dim = vFun.getDim();
		MathFunc rlt = C0;
		for(int i=1; i<=dim; i++) {
			MathFunc fd = (MathFunc)vFun.get(i);
			rlt = rlt.A(fd.diff(vFun.varNames().get(i-1)));
		}
		return rlt;
	}

	/**
	 * Compute divergence of <code>vFun</code> with respect to variables <code>varNames</code>
	 * in case of composition of functions.
	 * 
	 * @param fun
	 * @param varNames
	 * @return
	 */
	public static MathFunc div(VectorMathFunc vFun,String ...varNames) {
		MathFunc rlt = new FC(0.0);
		for(int i=0;i<varNames.length;i++) {
			MathFunc fd = (MathFunc)vFun.get(i+1);
			rlt = rlt.A(fd.diff(varNames[i]));
		}
		return rlt;
	}
	
	/**
	 * Compute divergence of <code>vFun</code> with respect to variables <code>varNames</code>
	 * in case of composition of functions.
	 * 
	 * @param fun
	 * @param varNames
	 * @return
	 */
	public static MathFunc div(VectorMathFunc vFun,ObjList<String> varNames) {
		MathFunc rlt = new FC(0.0);
		for(int i=1;i<=varNames.size();i++) {
			MathFunc fd = (MathFunc)vFun.get(i);
			rlt = rlt.A(fd.diff(varNames.at(i)));
		}
		return rlt;
	}
	
	public static MathFunc curl(VectorMathFunc vFun) {
		//TODO
		return null;
	}
	public static MathFunc curl(VectorMathFunc vFun, String ...varNames) {
		//TODO
		return null;
	}	
	public static MathFunc curl(VectorMathFunc vFun, ObjList<String> varNames) {
		//TODO
		return null;
	}
	
	
	//--- Vectors operations----------------------------------
	
	/**
	 * Vector summation
	 * 
	 * @param vi Vectors
	 * @return Vector result = v1 + v2 + ... +vn
	 */
	public static Vector sum(Vector ...vi) {
		if(vi==null || vi.length==0) {
			throw new FutureyeException("Check parameters vi="+vi);
		}
		Vector rlt = vi[0].copy();
		for(int i=1;i<vi.length;i++) {
			rlt = rlt.add(vi[i]);
		}
		return rlt;
	}
	
	/**
	 * Sum all the components of a vector
	 * @param v
	 * @return v_1 + v_2 + ... + v_n
	 */
	public static double sum(Vector v) {
		if(v == null || v.getDim() == 0)
			throw new FutureyeException("It should be at least one value in vector v!");
		double rlt = 0.0;
		int dim = v.getDim();
		for(int i=1;i<=dim;i++) {
			rlt += v.get(i);
		}
		return rlt;
	}
	
	/**
	 * Returns the natural logarithm (base e) of each component of vector v
	 * @param v
	 * @return
	 */
	public static Vector log(Vector v) {
		Vector v2 = v.copy();
		for(int i=1;i<=v.getDim();i++) {
			v2.set(i,Math.log(v.get(i)));
		}
		return v2;
	}
	
	/**
	 * Returns the base 10 logarithm of each component of vector v
	 * @param v
	 * @return
	 */
	public static Vector log10(Vector v) {
		Vector v2 = v.copy();
		for(int i=1;i<=v.getDim();i++) {
			v2.set(i,Math.log10(v.get(i)));
		}
		return v2;		
	}
	
	public static Vector exp(Vector v) {
		Vector v2 = v.copy();
		for(int i=1;i<=v.getDim();i++) {
			v2.set(i,Math.exp(v.get(i)));
		}
		return v2;
	}
	
	public static Vector abs(Vector v) {
		Vector v2 = v.copy();
		for(int i=1;i<=v.getDim();i++) {
			v2.set(i,Math.abs(v.get(i)));
		}
		return v2;
	}
	
	public static double max(Vector v) {
		if(v == null || v.getDim() == 0)
			throw new FutureyeException("It should be at least one value in vector v!");
		double max = v.get(1);
		for(int i=2;i<=v.getDim();i++) {
			double val = v.get(i);
			if(val > max) max = val;
		}
		return max;
	}
	
	public static double max(double[] a) {
		if(a == null || a.length == 0)
			throw new FutureyeException("It should be at least one value in array a!");
		double max = a[0];
		for(int i=1;i<a.length;i++) {
			if(a[i] > max) max = a[i];
		}
		return max;
	}
	
	public static double min(Vector v) {
		if(v == null || v.getDim() == 0)
			throw new FutureyeException("It should be at least one value in vector v!");
		double min = v.get(1);
		for(int i=2;i<=v.getDim();i++) {
			double val = v.get(i);
			if(val < min) min = val; 
		}
		return min;
	}
	
	public static double min(double[] a) {
		if(a == null || a.length == 0)
			throw new FutureyeException("It should be at least one value in array a!");
		double min = a[0];
		for(int i=0;i<a.length;i++) {
			if(a[i] < min) min = a[i];
		}
		return min;
	}
	
	/**
	 * y=a*x
	 * @param a
	 * @param x
	 * @return
	 */
	public static Vector ax(double a, Vector x) {
		Vector rlt = x.copy();
		return rlt.ax(a);
	}
	
	/**
	 * z = a*x+y
	 * @param a
	 * @param x
	 * @param y
	 * @return
	 */
	public static Vector axpy(double a, Vector x, Vector y) {
		Vector rlt = x.copy();
		return rlt.axpy(a, y);
	}
	
	/**
	 * zi = a*xi*yi 
	 * @param a
	 * @param x
	 * @param y
	 * @return
	 */
	public static Vector axMuly(double a, Vector x, Vector y) {
		Vector rlt = x.copy();
		return rlt.axMuly(a, y);
	}
	
	/**
	 * zi = a*xi/yi
	 * @param a
	 * @param x
	 * @param y
	 * @return
	 */
	public static Vector axDivy(double a, Vector x, Vector y) {
		Vector rlt = x.copy();
		return rlt.axDivy(a, y);

	}
	
	/**
	 * xi^exp
	 * 
	 * @param x
	 * @param exp
	 * @return
	 */
	public static Vector pow(Vector x, double exp) {
		Vector rlt = x.copy();
		for(int i=1;i<=x.getDim();i++) {
			rlt.set(i,Math.pow(x.get(i), exp));
		}
		return rlt;
	}
	
	/**
	 * base^xi
	 * 
	 * @param base
	 * @param x
	 * @return
	 */
	public static Vector pow(double base, Vector x) {
		Vector rlt = x.copy();
		for(int i=1;i<=x.getDim();i++) {
			rlt.set(i,Math.pow(base, x.get(i)));
		}
		return rlt;
		
	}
	
	//------statistic functions--------------------------------
	
	public static double mean(Vector v) {
		double mn = 0.0;
		for(int i=v.getDim()+1; --i>0;)
			mn += v.get(i);
		mn /= v.getDim();
		return mn;
	}
	
	public static double variance(Vector v) {
		double mnv = mean(v);
		double varv = 0.0, tmp;
		for(int i=v.getDim()+1; --i>0;) {
			tmp = mnv - v.get(i);
			varv += tmp*tmp;
		}
		varv /= v.getDim();
		return varv;
	}
	
	/**
	 * Standard Deviation
	 * @param v
	 * @return
	 */
	public static double SD(Vector v) {
		return Math.sqrt(variance(v));
	}
	
	/**
	 * Sample Standard Deviation
	 * @param v
	 * @return
	 */
	public static double sampleSD(Vector v) {
		double mnv = mean(v);
		double varv = 0.0, tmp;
		for(int i=v.getDim()+1; --i>0;) {
			tmp = mnv - v.get(i);
			varv += tmp*tmp;
		}
		varv /= v.getDim()-1;
		return varv;
	}
	
	public static double averageAbsoluteDeviation(Vector v) {
		double mnv = mean(v);
		double aad = 0.0;
		for(int i=v.getDim()+1; --i>0;) {
			aad += Math.abs(mnv - v.get(i));
		}
		aad /= v.getDim();
		return aad;
	}

}
