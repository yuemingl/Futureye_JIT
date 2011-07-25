package edu.uta.futureye.function.operator;

import java.util.List;

import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.function.AbstractFunction;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.basic.SpaceVectorFunction;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.function.intf.VectorFunction;
import edu.uta.futureye.util.FutureyeException;
import edu.uta.futureye.util.Utils;
import edu.uta.futureye.util.container.ObjList;

public class FMath {
	
	////////////////////////////Basic////////////////////////
	public static Function sqrt(final Function f) {
		return new AbstractFunction(f.varNames()) {
			@Override
			public double value(Variable v) {
				return Math.sqrt(f.value(v));
			}
			@Override
			public Function _d(String varName) {
				return FC.c(0.5).M(pow(f,-0.5)).M(f._d(varName));
			}
			@Override
			public int getOpOrder() {
				return OP_ORDER1;
			}
			@Override
			public String toString() {
				return "sqrt("+f.toString()+")";
			}
		};
	}
	
	public static Function pow(final Function f, final double p) {
		return new AbstractFunction(f.varNames()) {
			@Override
			public double value(Variable v) {
				return Math.pow(f.value(v),p);
			}
			@Override
			public Function _d(String varName) {
				return FC.c(p).M(pow(f,p-1)).M(f._d(varName));
			}
			@Override
			public int getOpOrder() {
				return OP_ORDER1;
			}
			@Override
			public String toString() {
				if( (f instanceof FC && Double.compare(f.value(), 0.0)==0))
					return "";
				return "("+f.toString()+")^"+p+"";
			}
		};
	}

	public static Function pow(final Function f1, final Function f2) {
		return new AbstractFunction(Utils.mergeList(f1.varNames(), f2.varNames())) {
			@Override
			public double value(Variable v) {
				return Math.pow(f1.value(v),f2.value(v));
			}
			@Override
			public int getOpOrder() {
				return OP_ORDER1;
			}
			@Override
			public String toString() {
				if( (f1 instanceof FC && Double.compare(f1.value(), 0.0)==0))
					return "";
				return "("+f1.toString()+")^("+f2.toString()+")";
			}
		};
	}	
	
	public static Function abs(final Function f) {
		return new AbstractFunction(f.varNames()) {
			@Override
			public double value(Variable v) {
				return Math.abs(f.value(v));
			}
			@Override
			public int getOpOrder() {
				return OP_ORDER1;
			}
			@Override
			public String toString() {
				return "abs("+f.toString()+")";
			}
		};
	}
	
	/**
	 * \sum{f_i} = f_1 + f_2 + ... + f_N
	 * 
	 * @param fi
	 * @return
	 */
	public static Function sum(Function ...fi) {
		if(fi==null || fi.length==0) {
			throw new FutureyeException("Check parameter f="+fi);
		}
		Function rlt = fi[0];
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
	public static Function linearCombination(double c1, Function f1,
			double c2,Function f2) {
		return new FC(c1).M(f1).A(new FC(c2).M(f2));
	}
	
	/**
	 * \sum{c_i*f_i} = c_1*f_1 + c_2*f_2 + ... + c_N*f_N
	 * 
	 * @param ci
	 * @param fi
	 * @return
	 */
	public static Function linearCombination(double []ci, Function []fi) {
		if(ci==null || ci.length==0 || fi==null || fi.length==0) {
			throw new FutureyeException("Check parameters ci="+ci+", fi="+fi);
		} else if(ci.length != fi.length) {
			throw new FutureyeException("(ci.length="+ci.length+") != (fi.lenght="+fi.length+")");
		}
		Function rlt = new FC(ci[0]).M(fi[0]);
		for(int i=1;i<fi.length;i++) {
			rlt = rlt.A(new FC(ci[i]).M(fi[i]));
		}
		return rlt;
	}
	
	public static VectorFunction grad(Function fun) {
		List<String> names = fun.varNames();
		VectorFunction rlt = new SpaceVectorFunction(names.size());
		for(int i=0;i<names.size();i++)
			rlt.set(i+1, fun._d(names.get(i)));
		return rlt;
	}
	
	public static VectorFunction grad(Function fun,ObjList<String> varNames) {
		VectorFunction rlt = new SpaceVectorFunction(varNames.size());
		for(int i=1;i<=varNames.size();i++)
			rlt.set(i, fun._d(varNames.at(i)));
		return rlt;
	}
	
	public static Function div(VectorFunction vFun) {
		int dim = vFun.getDim();
		Function rlt = FC.c0;
		for(int i=1; i<=dim; i++) {
			Function fd = (Function)vFun.get(i);
			rlt = rlt.A(fd._d(vFun.varNames().get(i-1)));
		}
		return rlt;
	}
	
	public static Function div(VectorFunction vFun,ObjList<String> varNames) {
		Function rlt = new FC(0.0);
		for(int i=1;i<=varNames.size();i++) {
			Function fd = (Function)vFun.get(i);
			rlt = rlt.A(fd._d(varNames.at(i)));
		}
		return rlt;
	}
	
	public static Function curl(VectorFunction vFun) {
		//TODO
		return null;
	}
	public static Function curl(VectorFunction vFun, ObjList<String> varNames) {
		//TODO
		return null;
	}
	
	/////////////////////Vectors///////////////////////
	
	public static Vector sum(Vector ...vi) {
		if(vi==null || vi.length==0) {
			throw new FutureyeException("Check parameter vi="+vi);
		}
		Vector rlt = vi[0].copy();
		for(int i=1;i<vi.length;i++) {
			rlt = rlt.add(vi[i]);
		}
		return rlt;
	}
	
	public static Vector log(Vector v) {
		Vector v2 = v.copy();
		for(int i=1;i<=v.getDim();i++) {
			v2.set(i,Math.log(v.get(i)));
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
		double max = Double.MIN_VALUE;
		for(int i=1;i<=v.getDim();i++) {
			double val = v.get(i);
			if(val > max)
				max = val; 
		}
		return max;
	}
	
	public static double min(Vector v) {
		double min = Double.MAX_VALUE;
		for(int i=1;i<=v.getDim();i++) {
			double val = v.get(i);
			if(val < min)
				min = val; 
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
}
