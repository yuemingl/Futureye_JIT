package edu.uta.futureye.function.operator;

import java.util.Map;
import java.util.Map.Entry;

import edu.uta.futureye.function.DerivativeIndicator;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.FAbstract;
import edu.uta.futureye.function.basic.FConstant;
import edu.uta.futureye.function.intf.FunctionDerivable;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.util.Utils;

public class FOBasicDerivable {
	public static FunctionDerivable Plus(final FunctionDerivable f1, final FunctionDerivable f2) {
		if(f1 instanceof FConstant && Double.compare(f1.value(null), 0.0)==0)
			return f2;
		else if(f2 instanceof FConstant && Double.compare(f2.value(null), 0.0)==0)
			return f1;
		else
			return new FAbstract(Utils.mergeList(f1.varNames(), f2.varNames())) {
				@Override
				public double value(Variable v) {
					return f1.value(v) + f2.value(v);
				}
				@Override
				public FunctionDerivable derivative(DerivativeIndicator di) {
					return Plus(f1.derivative(di),f2.derivative(di));
				}
				public String toString() {
					return "("+f1.toString()+"  +  "+f2.toString()+")";
				}
			};
	}
	
	public static FunctionDerivable Minus(final FunctionDerivable f1, final FunctionDerivable f2) {
		if(f1 instanceof FConstant && Double.compare(f1.value(null), 0.0)==0 &&
				f2 instanceof FConstant)
			return new FConstant( - f2.value(null) );
		else if(f2 instanceof FConstant && Double.compare(f2.value(null), 0.0)==0)
			return f1;
		else
			return new FAbstract(Utils.mergeList(f1.varNames(), f2.varNames())) {
				@Override
				public double value(Variable v) {
					return f1.value(v) - f2.value(v);
				}
				@Override
				public FunctionDerivable derivative(DerivativeIndicator di) {
					return Minus(f1.derivative(di),f2.derivative(di));
				}
				public String toString() {
					if(f1 instanceof FConstant && Double.compare(f1.value(null), 0.0)==0)
						return "  -  "+f2.toString();
					return "("+f1.toString()+"  -  "+f2.toString()+")";
				}
			};
	}
	
	public static FunctionDerivable Mult(final FunctionDerivable f1, final FunctionDerivable f2) {
		if( (f1 instanceof FConstant && Double.compare(f1.value(null), 0.0)==0) ||
				f2 instanceof FConstant && Double.compare(f2.value(null), 0.0)==0)
			return new FConstant(0.0);
		else if(f1 instanceof FConstant && Double.compare(f1.value(null), 1.0)==0)
			return f2;
		else if(f2 instanceof FConstant && Double.compare(f2.value(null), 1.0)==0)
			return f1;
		else
			return new FAbstract(Utils.mergeList(f1.varNames(), f2.varNames())) {
				@Override
				public double value(Variable v) {
					return f1.value(v) * f2.value(v);
				}
				@Override
				public FunctionDerivable derivative(DerivativeIndicator di) {
					return Plus(
							Mult(f1.derivative(di),f2),
							Mult(f1,f2.derivative(di))
							);
				}
				public String toString() {
					return "("+f1.toString()+"  <*>  "+f2.toString()+")";
				}
			};
	}
	
	public static FunctionDerivable Divi(final FunctionDerivable f1, final FunctionDerivable f2) {
		if( (f1 instanceof FConstant && Double.compare(f1.value(null), 0.0)==0))
			return new FConstant(0.0);
		return new FAbstract(Utils.mergeList(f1.varNames(), f2.varNames())) {
			@Override
			public double value(Variable v) {
				return f1.value(v) / f2.value(v);
			}
			@Override
			public FunctionDerivable derivative(DerivativeIndicator di) {
				return Plus(
						Mult(f1.derivative(di),f2),
						Mult(f1,f2.derivative(di))
						);
			}
			public String toString() {
				return "("+f1.toString()+"  </>  "+f2.toString()+")";
			}
		};
	}
	
	public static FunctionDerivable LinearCombination(final double coef1, final FunctionDerivable f1,
			final double coef2,final FunctionDerivable f2) {
		return new FAbstract(Utils.mergeList(f1.varNames(), f2.varNames())) {
			@Override
			public double value(Variable v) {
				return coef1*f1.value(v) + coef2 * f2.value(v);
			}
			@Override
			public FunctionDerivable derivative(DerivativeIndicator di) {
				return Plus(
						Mult(new FConstant(coef1),f1.derivative(di)),
						Mult(new FConstant(coef2),f2.derivative(di))
						);
			}
			public String toString() {
				return "("+coef1+"*("+f1.toString()+")  +  "+coef2+"*("+f2.toString()+"))";
			}
		};
	}
	
	/**
	 * 复合函数
	 * @param e.g. fOuter f = f(x,y)
	 * @param e.g. fInners x = x(r,s), y = y(r,s)
	 * @return e.g. f = f(x,y) = f( x(r,s),y(r,s) )
	 */
	public static FunctionDerivable Compose(final FunctionDerivable fOuter, final Map<String,FunctionDerivable> fInners) {
		return new FAbstract(fOuter.varNames()) {
			@Override
			public double value(Variable v) {
				//bugfix 增加或条件
				if(fOuter.varNames().containsAll(v.getValues().keySet()) ||
						v.getValues().keySet().containsAll(fOuter.varNames())) {
					return fOuter.value(v);
				} else if(fOuter.varNames().size() == fInners.size()){
					Variable newVar = new Variable();
					for(String varName : fOuter.varNames()) {
						Function f = fInners.get(varName);
						if(f != null ) {
							newVar.set(varName, f.value(v));
						}
						else {
							System.out.println("ERROR: Can not find "+varName+" in fInners.");
							return 0.0;
						}
					}
					return fOuter.value(newVar);
				} else {
					System.out.println("ERROR: Variable Number Mismatch of fOuter and fInner.");
				}
				return 0.0;
			}
			/**
			 * 链式求导
			 * f( x(r,s),y(r,s) )_r = f_x * x_r + f_y * y_r
			 */
			@Override
			public FunctionDerivable derivative(DerivativeIndicator di) {
				Map<String,Integer> diMap = di.get();
				for(Entry<String, Integer> diEntry : diMap.entrySet()) {
					
					if(fOuter.varNames().contains(diEntry.getKey())) {
						//f(x,y)关于x或y求导
						FunctionDerivable rlt = fOuter.derivative(di);
						return rlt;
					} else {
						//f(x,y)关于r或s求导
						FunctionDerivable rlt = new FConstant(0.0);
						for(String varName : fOuter.varNames()) {
							DerivativeIndicator diOut = new DerivativeIndicator();
							//TODO
							diOut.set(varName, 1);
							FunctionDerivable fInner = fInners.get(varName);
							if(fInner != null) {
								FunctionDerivable rltOuter = fOuter.derivative(diOut);
								FunctionDerivable rltInner = fInner.derivative(di);
								//f_x * x_r + f_y * y_r
								rlt = Plus(rlt, Mult(rltOuter,rltInner));
							}
						}
						return rlt;
					}
				}
				return null;
				
			}
			public String toString() {
				return "("+fOuter.toString()+"  ,  "+fInners.toString()+")";
			}
		};
	}
	
	public static FunctionDerivable Sqrt(final Function f) {
		return new FAbstract(f.varNames()) {
			@Override
			public double value(Variable v) {
				return Math.sqrt(f.value(v));
			}
			@Override
			public FunctionDerivable derivative(DerivativeIndicator di) {
				return null;
			}
			public String toString() {
				return "Sqrt("+f.toString()+")";
			}
			
		};
	}
	
	public static FunctionDerivable Abs(final Function f) {
		return new FAbstract(f.varNames()) {
			@Override
			public double value(Variable v) {
				return Math.abs(f.value(v));
			}
			@Override
			public FunctionDerivable derivative(DerivativeIndicator di) {
				return null;
			}
			public String toString() {
				return "Abs("+f.toString()+")";
			}
		};
	}		
}
