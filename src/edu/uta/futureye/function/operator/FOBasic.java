package edu.uta.futureye.function.operator;

import java.util.List;

import edu.uta.futureye.function.AbstractFunction;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.util.FutureyeException;
import edu.uta.futureye.util.Utils;

public class FOBasic {
	public static Function Plus(final Function f1, final Function f2) {
		return new AbstractFunction(Utils.mergeList(f1.varNames(), f2.varNames())) {
			@Override
			public double value(Variable v) {
				return f1.value(v) + f2.value(v);
			}
			public String toString() {
				if(f1 instanceof FC && Double.compare(f1.value(null), 0.0)==0)
					return f2.toString();
				else if(f2 instanceof FC && Double.compare(f2.value(null), 0.0)==0)
					return f1.toString();
				return "("+f1.toString()+"  +  "+f2.toString()+")";
			}
		};
	}
	
	public static Function PlusAll(final Function ...f) {
		if(f.length<=1) {
			FutureyeException e = new FutureyeException("Parameter length should > 1");
			e.printStackTrace();
			return null;
		}
		List<String> names = f[0].varNames();
		for(int i=1;i<f.length;i++) {
			names = Utils.mergeList(names, f[i].varNames());
		}
		return new AbstractFunction(names) {
			@Override
			public double value(Variable v) {
				double val = f[0].value(v);
				for(int i=1;i<f.length;i++) 
					val += f[i].value(v);
				return val;
			}
			public String toString() {
				String name = "("+f[0].toString();
				for(int i=1;i<f.length;i++) 
					name += " + "+f[i].toString();
				name += ")";
				return name;
			}
		};
	}
	
	public static Function Minus(final Function f1, final Function f2) {
		return new AbstractFunction(Utils.mergeList(f1.varNames(), f2.varNames())) {
			@Override
			public double value(Variable v) {
				return f1.value(v) - f2.value(v);
			}
			public String toString() {
				if(f1 instanceof FC && Double.compare(f1.value(null), 0.0)==0)
					return "  -  "+f2.toString();
				else if(f2 instanceof FC && Double.compare(f2.value(null), 0.0)==0)
					return f1.toString();
				return "("+f1.toString()+"  -  "+f2.toString()+")";
			}
		};
	}
	
	public static Function Mult(final Function f1, final Function f2) {
		return new AbstractFunction(Utils.mergeList(f1.varNames(), f2.varNames())) {
			@Override
			public double value(Variable v) {
				return f1.value(v) * f2.value(v);
			}
			public String toString() {
				if( (f1 instanceof FC && Double.compare(f1.value(null), 0.0)==0) ||
						f2 instanceof FC && Double.compare(f2.value(null), 0.0)==0)
					return "";
				else if(f1 instanceof FC && Double.compare(f1.value(null), 1.0)==0)
					return f2.toString();
				else if(f2 instanceof FC && Double.compare(f2.value(null), 1.0)==0)
					return f1.toString();
				return "("+f1.toString()+"  <*>  "+f2.toString()+")";
			}	
		};
	}
	
	public static Function Divi(final Function f1, final Function f2) {
		return new AbstractFunction(Utils.mergeList(f1.varNames(), f2.varNames())) {
			@Override
			public double value(Variable v) {
				return f1.value(v) / f2.value(v);
			}
			public String toString() {
				if( (f1 instanceof FC && Double.compare(f1.value(null), 0.0)==0))
					return "";		
				return "("+f1.toString()+"  </>  "+f2.toString()+")";
			}			
		};
	}

	public static Function Sqrt(final Function f) {
		return new AbstractFunction(f.varNames()) {
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
		return new AbstractFunction(f.varNames()) {
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
		return new AbstractFunction(Utils.mergeList(f1.varNames(), f2.varNames())) {
			@Override
			public double value(Variable v) {
				return Math.pow(f1.value(v),f2.value(v));
			}
			public String toString() {
				if( (f1 instanceof FC && Double.compare(f1.value(null), 0.0)==0))
					return "";		
				return "("+f1.toString()+")^("+f2.toString()+")";
			}			
		};
	}	
}
