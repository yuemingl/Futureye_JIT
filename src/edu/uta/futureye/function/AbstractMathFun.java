package edu.uta.futureye.function;

import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.intf.MathFun;
import edu.uta.futureye.function.intf.SimpleFun;
import edu.uta.futureye.util.Constant;
import edu.uta.futureye.util.FutureyeException;
import edu.uta.futureye.util.Utils;

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
		final MathFun fOuter = this;
		boolean find = false;
		for(String key : fInners.keySet()) {
			if(getVarNames().contains(key)) find = true;
		}
		if(!find) 
			return this; //No compose
		else
			return new AbstractMathFun(fOuter.getVarNames()) {
			@Override
			public double apply(Variable v) {
				return apply(v,null);
			}
			
			@Override
			public double apply(Variable v, Map<Object,Object> cache) {
				
				//if(fOuter.varNames().size() == 0) {
				//	throw new FutureyeException("\nERROR:\n fOuter varNames list is empty!");
				//}
				
				//bugfix 增加或条件 
				//bugfix 3/19/12
				//bug?3/20/12  v=[r], fOuter.varNames()=[s,t], 但fOuter的表达式只有r, 这种情况下会进入else分支，
				//一般来说是不会有这种情况的，如果确实有这种情况，需要在函数类增加activeVarNames
				//if(fOuter.varNames().containsAll(v.getValues().keySet()) ||
				//		v.getValues().keySet().containsAll(fOuter.varNames())) {
				if(v.getNameValuePairs().keySet().containsAll(fOuter.getVarNames())) {
					return fOuter.apply(v,cache);
				//} else if(fOuter.varNames().size() == fInners.size()){
				} else {
					Variable newVar = new Variable();
					for(String varName : fOuter.getVarNames()) {
						MathFun fInner = fInners.get(varName);
						if(fInner != null ) 
							newVar.set(varName, fInner.apply(v,cache));
						else //for mixed case: fOuter( x(r,s,t), y(r,s,t), r, s) bugfix 3/19/12
							newVar.set(varName, v.get(varName));
						//	throw new FutureyeException("\nERROR:\n Can not find "+varName+" in fInners.");
					}
					return fOuter.apply(newVar,cache);
				}
//				else {
//					throw new FutureyeException(
//							"\nERROR:\n Variable number mismatch of fOuter("+
//							fOuter.varNames()+") and fInner("+fInners+").");
//				}
			} 

			
			@Override
			public double[] applyAll(VariableArray v, Map<Object,Object> cache) {
				//bugfix 增加或条件
				if(v.getValues().keySet().containsAll(fOuter.getVarNames())) {
					return fOuter.applyAll(v,cache);
				} else {
					VariableArray newVar = new VariableArray();
					for(String varName : fOuter.getVarNames()) {
						MathFun fInner = fInners.get(varName);
						if(fInner != null )
							newVar.set(varName, fInner.applyAll(v,cache));
						else //for mixed case: fOuter( x(r,s,t), y(r,s,t), r, s)
							newVar.set(varName, v.get(varName));
					}
					return fOuter.applyAll(newVar,cache);
				}
			}
			
			/**
			 * 链式求导
			 * f( x(r,s),y(r,s) )_r = f_x * x_r + f_y * y_r
			 */
			@Override
			public MathFun _d(String varName) {
				MathFun rlt = null;
				if(fOuter.getVarNames().contains(varName)) {
					//f(x,y)关于x或y求导
					rlt = fOuter._d(varName);
					return rlt;
				} else {
					//f(x,y)关于r或s求导
					rlt = new FC(0.0);
					for(String innerVarName : fOuter.getVarNames()) {
						MathFun fInner = fInners.get(innerVarName);
						if(fInner != null) {
							MathFun rltOuter = fOuter._d(innerVarName);
							if(!(rltOuter.isConstant()))
								rltOuter = rltOuter.compose(fInners);
							MathFun rltInner = fInner._d(varName);
							//f_x * x_r + f_y * y_r
							rlt = rlt.A(
									rltOuter.M(rltInner)
									);
						}
					}
					return rlt;
				}
			}
			@Override
			public int getOpOrder() {
				return fOuter.getOpOrder();
			}
			@Override
			public String toString() {
				String rlt = fOuter.toString();
				for(Entry<String,MathFun> map : fInners.entrySet()) {
					String names = map.getValue().getVarNames().toString();
					rlt = rlt.replace(map.getKey(), 
							map.getKey()+"("+names.substring(1,names.length()-1)+")");
				}
				return rlt;
			}
		};
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
			return new AbstractMathFun(Utils.mergeList(f1.getVarNames(), f2.getVarNames())) {
				@Override
				public double apply(Variable v) {
					return f1.apply(v) + f2.apply(v);
				}
				
				@Override
				public double apply(Variable v, Map<Object,Object> cache) {
//基本运算不需要cache，否则计算效率会更低					
//					if(cache != null) {
//						Double v1, v2;
//						v1 = cache.get(f1);
//						if(v1 == null) {
//							v1 = f1.value(v,cache);
//							cache.put(f1, v1);
//						}
//						v2 = cache.get(f2);
//						if(v2 == null) {
//							v2 = f2.value(v,cache);
//							cache.put(f2, v2);
//						}
//						return v1 + v2;
//					} else {
//						return value(v);
//					}
					return f1.apply(v,cache) + f2.apply(v,cache);

				}
				
				@Override
				public double[] applyAll(VariableArray v, Map<Object,Object> cache) {
					int len = v.length();
					double[] la = f1.applyAll(v,cache);
					double[] ra = f2.applyAll(v,cache);
					for(int i=0;i<len;i++) {
						la[i] += ra[i];
					}
					return la;
				}
				
				@Override
				public MathFun _d(String varName) {
					return f1._d(varName).A(f2._d(varName)).setVarNames(this.getVarNames());
				}
				@Override
				public int getOpOrder() {
					return OP_ORDER3;
				}
				@Override
				public String toString() {
					StringBuilder sb = new StringBuilder();
					sb.append(f1.toString());
					sb.append(" + ");
					sb.append(f2.toString());
					return sb.toString();
				}
			};
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
			return new AbstractMathFun(Utils.mergeList(f1.getVarNames(), f2.getVarNames())) {
				@Override
				public double apply(Variable v) {
					return f1.apply(v) - f2.apply(v);
				}
				
				@Override
				public double apply(Variable v, Map<Object,Object> cache) {
					return f1.apply(v,cache) - f2.apply(v,cache);
				}
				
				@Override
				public double[] applyAll(VariableArray v, Map<Object,Object> cache) {
					int len = v.length();
					double[] la = f1.applyAll(v,cache);
					double[] ra = f2.applyAll(v,cache);
					for(int i=0;i<len;i++) {
						la[i] -= ra[i];
					}
					return la;
				}
				
				@Override
				public MathFun _d(String varName) {
					return f1._d(varName).S(f2._d(varName)).setVarNames(this.getVarNames());
				}
				@Override
				public int getOpOrder() {
					return OP_ORDER3;
				}
				@Override
				public String toString() {
					StringBuilder sb = new StringBuilder();
					if(! (f1.isConstant() && Math.abs(f1.apply()) < Constant.eps) ) {
						sb.append(f1.toString());
					}
					sb.append(" - ");
					if(f2.getOpOrder() >= OP_ORDER3)
						sb.append("(").append(f2.toString()).append(")");
					else
						sb.append(f2.toString());
					return sb.toString();
				}
			};
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
			return new AbstractMathFun(Utils.mergeList(f1.getVarNames(), f2.getVarNames())) {
				@Override
				public double apply(Variable v) {
					return f1.apply(v) * f2.apply(v);
				}
				
				@Override
				public double apply(Variable v, Map<Object,Object> cache) {
					return f1.apply(v,cache) * f2.apply(v,cache);

				}
				
				@Override
				public double[] applyAll(VariableArray v, Map<Object,Object> cache) {
					int len = v.length();
					double[] la = f1.applyAll(v,cache);
					double[] ra = f2.applyAll(v,cache);
					for(int i=0;i<len;i++) {
						la[i] *= ra[i];
					}
					return la;
				}
				
				@Override
				public MathFun _d(String varName) {
					return 	f1._d(varName).M(f2).A(
							f1.M(f2._d(varName))
							).setVarNames(this.getVarNames());
				}
				@Override
				public int getOpOrder() {
					return OP_ORDER2;
				}
				@Override
				public String toString() {
					StringBuilder sb = new StringBuilder();
					if(f1.getOpOrder() > OP_ORDER2)
						sb.append("(").append(f1.toString()).append(")");
					else
						sb.append(f1.toString());
					sb.append(" * ");
					if(f2.getOpOrder() > OP_ORDER2)
						sb.append("(").append(f2.toString()).append(")");
					else
						sb.append(f2.toString());
					return sb.toString();
				}
			};
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
			return new AbstractMathFun(Utils.mergeList(f1.getVarNames(), f2.getVarNames())) {
				@Override
				public double apply(Variable v) {
					return f1.apply(v) / f2.apply(v);
				}
				
				@Override
				public double apply(Variable v, Map<Object,Object> cache) {
					return f1.apply(v,cache) / f2.apply(v,cache);
				}
				
				@Override
				public double[] applyAll(VariableArray v, Map<Object,Object> cache) {
					int len = v.length();
					double[] la = f1.applyAll(v,cache);
					double[] ra = f2.applyAll(v,cache);
					for(int i=0;i<len;i++) {
						la[i] /= ra[i];
					}
					return la;
				}
				
				@Override
				public MathFun _d(String varName) {
					return f1._d(varName).M(f2).S(f1.M(f2._d(varName)))
							.D(f2.M(f2)).setVarNames(this.getVarNames());
				}
				@Override
				public int getOpOrder() {
					return OP_ORDER2;
				}
				@Override
				public String toString() {
					StringBuilder sb = new StringBuilder();
					if(f1.getOpOrder() > OP_ORDER2)
						sb.append("(").append(f1.toString()).append(")");
					else
						sb.append(f1.toString());
					sb.append(" / ");
					if(f2.getOpOrder() >= OP_ORDER2) //!!!
						sb.append("(").append(f2.toString()).append(")");
					else
						sb.append(f2.toString());
					return sb.toString();
				}
			};
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
