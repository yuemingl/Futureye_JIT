package edu.uta.futureye.function;

import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.util.Constant;
import edu.uta.futureye.util.FutureyeException;
import edu.uta.futureye.util.Utils;

public abstract class AbstractFunction implements Function {
	protected List<String> varNames = new LinkedList<String>();

	/**
	 * Construct n Dim function
	 * 
	 * @param varName
	 * @param aryVarNames
	 */
	public AbstractFunction(String varName, String ...aryVarNames) {
		varNames.add(varName);
		for(String s : aryVarNames)
			varNames.add(s);
	}
	
	/**
	 * Construct n Dim function
	 * 
	 * @param varNames
	 */
	public AbstractFunction(List<String> varNames) {
		this.varNames = varNames;
	}
	
	/**
	 * Construct 1 Dim function
	 */
	public AbstractFunction() {
	}
	
	@Override
	public List<String> varNames() {
		return varNames;
	}

	@Override
	public Function setVarNames(List<String> varNames) {
		this.varNames = varNames;
		return this;
	}
	
	/**
	 * Implement this function yourself
	 */
	@Override
	public double value(Variable v) {
		throw new UnsupportedOperationException();
	}
	
	@Override
	public double value() {
		throw new UnsupportedOperationException();
	}
	
	/**
	 * Implement this function yourself if necessary
	 */
	@Override
	public Function _d(String varName) {
		throw new UnsupportedOperationException();
	}

	@Override
	public Function compose(final Map<String,Function> fInners) {
		final Function fOuter = this;
		return new AbstractFunction(fOuter.varNames()) {
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
						if(f != null )
							newVar.set(varName, f.value(v));
						else
							throw new FutureyeException("\nERROR:\n Can not find "+
									varName+" in fInners.");
					}
					return fOuter.value(newVar);
				} else {
					throw new FutureyeException(
							"\nERROR:\n Variable number mismatch of fOuter("+
							fOuter.varNames()+") and fInner("+fInners+").");
				}
			}
			/**
			 * 链式求导
			 * f( x(r,s),y(r,s) )_r = f_x * x_r + f_y * y_r
			 */
			@Override
			public Function _d(String varName) {
				Function rlt = null;
				if(fOuter.varNames().contains(varName)) {
					//f(x,y)关于x或y求导
					rlt = fOuter._d(varName);
					return rlt;
				} else {
					//f(x,y)关于r或s求导
					rlt = new FC(0.0);
					for(String innerVarName : fOuter.varNames()) {
						Function fInner = fInners.get(innerVarName);
						if(fInner != null) {
							Function rltOuter = fOuter._d(innerVarName);
							if(!(rltOuter instanceof FC))
								rltOuter = rltOuter.compose(fInners);
							Function rltInner = fInner._d(varName);
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
				for(Entry<String,Function> map : fInners.entrySet()) {
					String names = map.getValue().varNames().toString();
					rlt = rlt.replace(map.getKey(), 
							map.getKey()+"("+names.substring(1,names.length()-1)+")");
				}
				return rlt;
			}
		};
	}

	
	////////////////////////Operations////////////////////////////////////
	
	@Override
	public Function A(final Function g) {
		final Function f1 = this;
		final Function f2 = g;
		if(f1 instanceof FC && f2 instanceof FC) {
			return new FC(f1.value() + f2.value());
		} else if(f1 instanceof FC && Math.abs(f1.value()) < Constant.eps) {
			return f2;
		} else if(f2 instanceof FC && Math.abs(f2.value()) < Constant.eps) {
			return f1;
		} else {
			return new AbstractFunction(Utils.mergeList(f1.varNames(), f2.varNames())) {
				@Override
				public double value(Variable v) {
					return f1.value(v) + f2.value(v);
				}
				@Override
				public Function _d(String varName) {
					return f1._d(varName).A(f2._d(varName)).setVarNames(this.varNames);
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
	public Function A(final double g) {
		return A(FC.c(g));
	}
	
	@Override
	public Function S(Function g) {
		final Function f1 = this;
		final Function f2 = g;
		if(f1 instanceof FC && f2 instanceof FC) {
			return new FC(f1.value() - f2.value());
		} else if(f2 instanceof FC && Math.abs(f2.value()) < Constant.eps) {
			return f1;
		} else {
			return new AbstractFunction(Utils.mergeList(f1.varNames(), f2.varNames())) {
				@Override
				public double value(Variable v) {
					return f1.value(v) - f2.value(v);
				}
				@Override
				public Function _d(String varName) {
					return f1._d(varName).S(f2._d(varName)).setVarNames(this.varNames);
				}
				@Override
				public int getOpOrder() {
					return OP_ORDER3;
				}
				@Override
				public String toString() {
					StringBuilder sb = new StringBuilder();
					if(! (f1 instanceof FC && Math.abs(f1.value()) < Constant.eps) ) {
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
	public Function S(final double g) {
		return S(FC.c(g));
	}	
	
	@Override
	public Function M(Function f) {
		final Function f1 = this;
		final Function f2 = f;
		if(f1 instanceof FC && f2 instanceof FC) {
			return new FC(f1.value() * f2.value());
		} else if( (f1 instanceof FC && Math.abs(f1.value()) < Constant.eps) ||
				f2 instanceof FC && Math.abs(f2.value()) < Constant.eps)
			return FC.c0;
		else if(f1 instanceof FC && Math.abs(f1.value()-1.0) < Constant.eps)
			return f2;
		else if(f2 instanceof FC && Math.abs(f2.value()-1.0) < Constant.eps)
			return f1;
		else
			return new AbstractFunction(Utils.mergeList(f1.varNames(), f2.varNames())) {
				@Override
				public double value(Variable v) {
					return f1.value(v) * f2.value(v);
				}
				@Override
				public Function _d(String varName) {
					return 	f1._d(varName).M(f2).A(
							f1.M(f2._d(varName))
							).setVarNames(this.varNames);
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
	public Function M(final double g) {
		return M(FC.c(g));
	}	
	
	@Override
	public Function D(Function f) {
		final Function f1 = this;
		final Function f2 = f;
		if(f1 instanceof FC && f2 instanceof FC) {
			return new FC(f1.value() / f2.value());
		} else if(f1 instanceof FC && Double.compare(f1.value(),0.0)==0) {
			//Math.abs(f1.value())<Constant.eps will not work properly
			return FC.c0;
		} else if(f2 instanceof FC && Double.compare(f2.value(),0.0)==0) {
			return FC.c(Double.POSITIVE_INFINITY);
		}  else if(f2 instanceof FC && Math.abs(f2.value()-1.0) < Constant.eps) {
			return f1;
		} else {
			return new AbstractFunction(Utils.mergeList(f1.varNames(), f2.varNames())) {
				@Override
				public double value(Variable v) {
					return f1.value(v) / f2.value(v);
				}
				@Override
				public Function _d(String varName) {
					return f1._d(varName).M(f2).S(f1.M(f2._d(varName)))
							.D(f2.M(f2)).setVarNames(this.varNames);
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
	public Function D(final double g) {
		return D(FC.c(g));
	}
	
	@Override
	public Function copy() {
		throw new UnsupportedOperationException();
	}
	////////////////////////For printing expression////////////////////////////
	
	@Override
	public String getFName() {
		return null;
	}
	
	@Override
	public void setFName(String name) {
		throw new UnsupportedOperationException();
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
	public String toString() {
		String s = varNames.toString();
		return "F("+s.substring(1, s.length()-1)+")";
	}
}
