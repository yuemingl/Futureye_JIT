package edu.uta.futureye.function;

import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.util.Constant;
import edu.uta.futureye.util.FutureyeException;
import edu.uta.futureye.util.Utils;

public abstract class AbstractFunction implements Function {
	protected List<String> varNames = new LinkedList<String>();;

	/**
	 * Construct 1 Dim function
	 */
	public AbstractFunction() {
		varNames.add(Constant.x);
	}
	
	/**
	 * Construct n Dim function
	 * @param varName
	 * @param aryVarNames
	 */
	public AbstractFunction(String varName, String ...aryVarNames) {
		varNames.add(varName);
		for(String s : aryVarNames)
			varNames.add(s);
	}
	
	public AbstractFunction(List<String> varNames) {
		this.varNames = varNames;
	}
	
	@Override
	public double value(Variable v) {
		throw new UnsupportedOperationException();
	}
	
	@Override
	public double value() {
		throw new UnsupportedOperationException();
	}
	
	@Override
	public Function d(String varName) {
		throw new UnsupportedOperationException();
	}
	
	@Override
	public List<String> varNames() {
		return varNames;
	}

	@Override
	public void setVarNames(List<String> varNames) {
		this.varNames = varNames;
	}
	
	@Override
	public Function P(final Function f) {
		final Function f1 = this;
		final Function f2 = f;
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
				public Function d(String varName) {
					return f1.d(varName).P(f2.d(varName));
				}
				public String toString() {
					return "("+f1.toString()+" [+] "+f2.toString()+")";
				}
			};
		}
	}
	
	@Override
	public Function M(Function f) {
		final Function f1 = this;
		final Function f2 = f;
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
				public Function d(String varName) {
					return f1.d(varName).M(f2.d(varName));
				}
				public String toString() {
					if(f1 instanceof FC && Math.abs(f1.value()) < Constant.eps)
						return " [-] "+f2.toString();
					return "("+f1.toString()+" [-] "+f2.toString()+")";
				}
			};
		}
	}
	
	@Override
	public Function X(Function f) {
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
				public Function d(String varName) {
					return 	f1.d(varName).X(f2).P(
							f1.X(f2.d(varName))
							);
				}
				public String toString() {
					return "("+f1.toString()+" [*] "+f2.toString()+")";
				}
			};
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
				public Function d(String varName) {
					return f1.d(varName).X(f2).M(
						   f1.X(f2.d(varName)).D(
								   f2.X(f2))
							);
				}
				public String toString() {
					return "("+f1.toString()+" [/] "+f2.toString()+")";
				}
			};
		}
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
						if(f != null ) {
							newVar.set(varName, f.value(v));
						}
						else {
							FutureyeException e = new FutureyeException("ERROR: Can not find "+varName+" in fInners.");
							e.printStackTrace();
							System.exit(0);
							return 0.0;
						}
					}
					return fOuter.value(newVar);
				} else {
					FutureyeException e = new FutureyeException(
							"ERROR: Variable Number Mismatch of fOuter("+
							fOuter.varNames()+") and fInner("+fInners+").");
					e.printStackTrace();
					System.exit(0);
				}
				return 0.0;
			}
			/**
			 * 链式求导
			 * f( x(r,s),y(r,s) )_r = f_x * x_r + f_y * y_r
			 */
			@Override
			public Function d(String varName) {
				Function rlt = null;
				if(fOuter.varNames().contains(varName)) {
					//f(x,y)关于x或y求导
					rlt = fOuter.d(varName);
					return rlt;
				} else {
					//f(x,y)关于r或s求导
					rlt = new FC(0.0);
					for(String innerVarName : fOuter.varNames()) {
						Function fInner = fInners.get(innerVarName);
						if(fInner != null) {
							Function rltOuter = fOuter.d(innerVarName);
							Function rltInner = fInner.d(varName);
							//f_x * x_r + f_y * y_r
							rlt = rlt.P(
									rltOuter.X(rltInner)
									);
						}
					}
					return rlt;
				}
			}
			public String toString() {
				//return "("+fOuter.toString()+"  ,  "+fInners.toString()+")";
				return "("+fOuter.toString()+")";
			}
		};
	}

	public String toString() {
		String s = varNames.toString();
		return "F("+s.substring(1, s.length()-1)+")";
	}
}
