package edu.uta.futureye.lib.shapefun;

import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import edu.uta.futureye.core.CoordinateTransform;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.function.AbstractFunction;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.FAxpb;
import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.basic.FX;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.function.intf.ScalarShapeFunction;
import edu.uta.futureye.function.operator.FMath;
import edu.uta.futureye.util.FutureyeException;
import edu.uta.futureye.util.container.ObjList;

/**
 * 三角形局部坐标，二次函数
 * 
 * 3
 * | \
 * |  \
 * 6   5
 * |    \
 * |     \
 * 1--4-- 2
 * 
 * N = N(r,s,t) = N( r(x,y), s(x,y), t(x,y) )
 * N1 = (2*r-1)*r
 * N2 = (2*s-1)*s
 * N3 = (2*t-1)*t
 * N4 = 4*r*s
 * N5 = 4*s*t
 * N6 = 4*r*t
 * 
 * @author liuyueming
 *
 */
public class SFQuadraticLocal2D extends AbstractFunction implements ScalarShapeFunction {
	private int funIndex;
	private Function funCompose = null;
	private Function funOuter = null;
	private List<String> varNames = new LinkedList<String>();
	private ObjList<String> innerVarNames = null;

	private Element e = null;

	public SFQuadraticLocal2D(int funID) {
		funIndex = funID - 1;
		if(funID<1 || funID>6) {
			FutureyeException ex = new FutureyeException("ERROR: funID should be 1~6.");
			ex.printStackTrace();
			System.exit(-1);
		}
		
		varNames.add("r");
		varNames.add("s");
		//varNames.add("t");
		innerVarNames = new ObjList<String>("x","y");
		
		//复合函数
		Map<String, Function> fInners = new HashMap<String, Function>(4);
		
		for(final String varName : varNames) {
			fInners.put(varName, new AbstractFunction(innerVarNames.toList()) {
				
				protected CoordinateTransform trans = new CoordinateTransform(2);
				
				public Function _d(String var) {
					//Coordinate transform and Jacbian on element e
					List<Function> funs = trans.getTransformFunction(
							//两种变换都可以，二次的要慢一些
							//trans.getTransformShapeFunctionByElement(e)
							trans.getTransformLinear2DShapeFunction(e)
								);
					trans.setTransformFunction(funs);
					
					Function fx = funs.get(0);
					Function fy = funs.get(1);
					
//					Function fr = new FX("r");
//					Function fs = new FX("s");
//					Map<String, Function> fInner_t = new HashMap<String, Function>(4);
//					fInner_t.put("t", FOBasicDerivable.Minus(
//							new FConstant(1.0), 
//							FOBasicDerivable.Plus(fr, fs)
//							));
//					Function fxCompose = FOBasicDerivable.Compose(fx, fInner_t);
//					Function fyCompose = FOBasicDerivable.Compose(fy, fInner_t);
//					
//					List<String> independentVarNames = new LinkedList<String>();
//					independentVarNames.add("r");
//					independentVarNames.add("s");
//					fxCompose.setVarNames(independentVarNames);
//					fyCompose.setVarNames(independentVarNames);

					Function x_r = fx._d("r");
					Function x_s = fx._d("s");
					Function y_r = fy._d("r");
					Function y_s = fy._d("s");
					
					//Function jac = (Function) trans.getJacobian2D();
					Function jac = (Function) e.getJacobin();
					
					if(varName.equals("r")) {
						if(var.equals("x"))
							return FMath.Divi(y_s, jac);
						if(var.equals("y"))
							return FMath.Minus(new FC(0.0),FMath.Divi(x_s, jac));
					} else if(varName.equals("s")) {
						if(var.equals("x"))
							return FMath.Minus(new FC(0.0),FMath.Divi(y_r, jac));
						if(var.equals("y"))
							return FMath.Divi(x_r, jac);
					} 
//					else if(varName.equals("t")) {
//						//t = 1 - r - s
//						//t_x = -r_x - s_x
//						//t_y = -r_y - s_y
//						if(di.get().get("x") != null)
//							return FOBasicDerivable.LinearCombination(
//									-0.1, FOBasicDerivable.Divi(y_s, jac), 
//									-0.1, FOBasicDerivable.Divi(y_r, jac));
//						if(di.get().get("y") != null)
//							return FOBasicDerivable.LinearCombination(
//									-0.1, FOBasicDerivable.Divi(x_s, jac), 
//									-0.1, FOBasicDerivable.Divi(x_r, jac));
//					}
					return null;
				}
			});
		}
		
		Function fr = new FX("r");
		Function fs = new FX("s");
		Function ft = //new FX("t");
		FMath.Minus(
				new FC(1.0), 
				FMath.Plus(fr, fs));
				
				
		Function f2rm1 = new FAxpb("r",2.0,-1.0);
		Function f2sm1 = new FAxpb("s",2.0,-1.0);
		Function f2tm1 = //new FAxpb("t",2.0,-1.0);
			FMath.Minus(
					new FC(1.0), 
					FMath.LinearCombination(2.0, fr, 2.0, fs));


//		Map<String, Function> fInner_t = new HashMap<String, Function>(4);
//		fInner_t.put("t", FOBasicDerivable.Minus(
//				new FConstant(1.0), 
//				FOBasicDerivable.Plus(fr, fs)
//				));
//		Function ftCompose = FOBasicDerivable.Compose(ft, fInner_t);
//		Function f2tm1Compose = FOBasicDerivable.Compose(f2tm1, fInner_t);
		Function ftCompose = ft;
		Function f2tm1Compose = f2tm1;
		
		if(funIndex == 0)
			funOuter = FMath.Mult(f2rm1, fr);
		else if(funIndex == 1)
			funOuter = FMath.Mult(f2sm1, fs);
		else if(funIndex == 2)
			funOuter = FMath.Mult(f2tm1Compose, ftCompose);
		else if(funIndex == 3)
			funOuter = FMath.Mult(new FC(4.0), FMath.Mult(fr, fs));
		else if(funIndex == 4)
			funOuter = FMath.Mult(new FC(4.0), FMath.Mult(fs, ftCompose));
		else if(funIndex == 5)
			funOuter = FMath.Mult(new FC(4.0), FMath.Mult(fr, ftCompose));
		
		//设置funOuter的独立自变量名称
//		List<String> independentVarNames = new LinkedList<String>();
//		independentVarNames.add("r");
//		independentVarNames.add("s");
//		funOuter.setVarNames(independentVarNames);
		funOuter.setVarNames(varNames);
		//使用复合函数构造形函数
		funCompose = FMath.Compose(funOuter, fInners);
		
	}

	@Override
	public void asignElement(Element e) {
		this.e = e;
	}

	//TODO ??? 应该采用一维二次型函数
	ScalarShapeFunction sf1d1 = new SFLinearLocal1D(1);
	ScalarShapeFunction sf1d2 = new SFLinearLocal1D(2);
	@Override
	public ScalarShapeFunction restrictTo(int funIndex) {
		if(funIndex == 1) return sf1d1;
		else return sf1d2;
	}

	@Override
	public Function _d(String varName) {
		return funCompose._d(varName);
	}

	@Override
	public void setVarNames(List<String> varNames) {
		this.varNames = varNames;
	}

	@Override
	public double value(Variable v) {
		return funCompose.value(v);
	}

	@Override
	public List<String> varNames() {
		return varNames;
	}

	public String toString() {
		return "N"+(funIndex+1)+"( r,s,t )="+funOuter.toString();
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		for(int i=1;i<=6;i++) {
			SFQuadraticLocal2D s = new SFQuadraticLocal2D(i);
			System.out.println(s);
			System.out.println(s._d("r"));
		}
	}

	@Override
	public ObjList<String> innerVarNames() {
		return innerVarNames;
	}

}
