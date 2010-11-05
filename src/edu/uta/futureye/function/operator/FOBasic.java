package edu.uta.futureye.function.operator;

import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.FAbstract;
import edu.uta.futureye.function.basic.FConstant;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.util.Utils;

public class FOBasic {
	public static Function Plus(final Function f1, final Function f2) {
		return new FAbstract(Utils.mergeList(f1.varNames(), f2.varNames())) {
			@Override
			public double value(Variable v) {
				return f1.value(v) + f2.value(v);
			}
			public String toString() {
				if(f1 instanceof FConstant && Double.compare(f1.value(null), 0.0)==0)
					return f2.toString();
				else if(f2 instanceof FConstant && Double.compare(f2.value(null), 0.0)==0)
					return f1.toString();
				return "("+f1.toString()+"  +  "+f2.toString()+")";
			}
		};
	}
	
	public static Function Minus(final Function f1, final Function f2) {
		return new FAbstract(Utils.mergeList(f1.varNames(), f2.varNames())) {
			@Override
			public double value(Variable v) {
				return f1.value(v) - f2.value(v);
			}
			public String toString() {
				if(f1 instanceof FConstant && Double.compare(f1.value(null), 0.0)==0)
					return "  -  "+f2.toString();
				else if(f2 instanceof FConstant && Double.compare(f2.value(null), 0.0)==0)
					return f1.toString();
				return "("+f1.toString()+"  -  "+f2.toString()+")";
			}
		};
	}
	
	public static Function Mult(final Function f1, final Function f2) {
		return new FAbstract(Utils.mergeList(f1.varNames(), f2.varNames())) {
			@Override
			public double value(Variable v) {
				return f1.value(v) * f2.value(v);
			}
			public String toString() {
				if( (f1 instanceof FConstant && Double.compare(f1.value(null), 0.0)==0) ||
						f2 instanceof FConstant && Double.compare(f2.value(null), 0.0)==0)
					return "";
				else if(f1 instanceof FConstant && Double.compare(f1.value(null), 1.0)==0)
					return f2.toString();
				else if(f2 instanceof FConstant && Double.compare(f1.value(null), 1.0)==0)
					return f1.toString();
				return "("+f1.toString()+"  <*>  "+f2.toString()+")";
			}	
		};
	}
	
	public static Function Divi(final Function f1, final Function f2) {
		return new FAbstract(Utils.mergeList(f1.varNames(), f2.varNames())) {
			@Override
			public double value(Variable v) {
				return f1.value(v) / f2.value(v);
			}
			public String toString() {
				if( (f1 instanceof FConstant && Double.compare(f1.value(null), 0.0)==0))
					return "";		
				return "("+f1.toString()+"  </>  "+f2.toString()+")";
			}			
		};
	}

	public static Function Sqrt(final Function f) {
		return new FAbstract(f.varNames()) {
			@Override
			public double value(Variable v) {
				return Math.sqrt(f.value(v));
			}
			public String toString() {
				return "Sqrt("+f.toString()+")";
			}
		};
	}
	
	public static Function Abs(final Function f) {
		return new FAbstract(f.varNames()) {
			@Override
			public double value(Variable v) {
				return Math.abs(f.value(v));
			}
			public String toString() {
				return "Abs("+f.toString()+")";
			}
		};
	}	
	
	public static Function Power(final Function f1, final Function f2) {
		return new FAbstract(Utils.mergeList(f1.varNames(), f2.varNames())) {
			@Override
			public double value(Variable v) {
				return Math.pow(f1.value(v),f2.value(v));
			}
			public String toString() {
				if( (f1 instanceof FConstant && Double.compare(f1.value(null), 0.0)==0))
					return "";		
				return "("+f1.toString()+")^("+f2.toString()+")";
			}			
		};
	}	
}
