package edu.uta.futureye.function.operator;

import java.util.List;
import java.util.Map;

import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.basic.SpaceVectorFunction;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.function.intf.VectorFunction;
import edu.uta.futureye.util.container.ObjList;

public class FOVector {
	
	/**
	 * 复合函数
	 */
	public static VectorFunction Compose(final VectorFunction fOuter, 
			final Map<String,Function> fInners) {
		
		int dim = fOuter.getDim();
		VectorFunction rlt = new SpaceVectorFunction(fOuter.getDim());
		for(int i=0;i<dim;i++)
			rlt.set(i+1, 
					FMath.Compose(
							(Function)fOuter.get(i+1), 
							fInners));
		return rlt;
	}	
	
	public static VectorFunction ScalarProduct(double coef, VectorFunction fun) {
		int dim = fun.getDim();
		SpaceVectorFunction svf = new SpaceVectorFunction(dim);
		for(int i=1;i<=dim;i++) {
			svf.set(i, FMath.Mult(
					new FC(coef),
					(Function)fun.get(i)));
		}
		return svf;
	}
	
	public static VectorFunction Plus(VectorFunction f1, VectorFunction f2) {
		int dim = f1.getDim();
		SpaceVectorFunction svf = new SpaceVectorFunction(dim);
		for(int i=1;i<=dim;i++) {
			svf.set(i, FMath.Plus(
					(Function)f1.get(i),
					(Function)f2.get(i)));
		}
		return svf;
	}
	
	public static VectorFunction Minus(VectorFunction f1, VectorFunction f2) {
		int dim = f1.getDim();
		SpaceVectorFunction svf = new SpaceVectorFunction(dim);
		for(int i=1;i<=dim;i++) {
			svf.set(i, FMath.Minus(
					(Function)f1.get(i),
					(Function)f2.get(i)));
		}
		return svf;
	}
	
	public static VectorFunction Grad(Function fun) {
		List<String> names = fun.varNames();
		VectorFunction rlt = new SpaceVectorFunction(names.size());
		for(int i=0;i<names.size();i++)
			rlt.set(i+1, fun._d(names.get(i)));
		return rlt;
	}
	
	public static VectorFunction Grad(Function fun,ObjList<String> varNames) {
		VectorFunction rlt = new SpaceVectorFunction(varNames.size());
		for(int i=1;i<=varNames.size();i++)
			rlt.set(i, fun._d(varNames.at(i)));
		return rlt;
	}
	
	public static Function Div(VectorFunction vFun) {
		int dim = vFun.getDim();
		Function rlt = new FC(0.0);
		for(int i=0; i<dim; i++) {
			Function fd = (Function)vFun.get(i+1);
			rlt = FMath.Plus(rlt, fd._d(vFun.varNames().get(i)));
		}
		return rlt;
	}
	
	public static Function Div(VectorFunction vFun,ObjList<String> varNames) {
		Function rlt = new FC(0.0);
		for(int i=1;i<=varNames.size();i++) {
			Function fd = (Function)vFun.get(i);
			rlt = FMath.Plus(rlt, fd._d(varNames.at(i)));
		}
		return rlt;
	}
	
	public static Function Curl(VectorFunction vFun) {
		//TODO
		return null;
	}
	
	
}
