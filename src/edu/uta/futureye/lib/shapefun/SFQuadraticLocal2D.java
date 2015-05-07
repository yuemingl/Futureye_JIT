package edu.uta.futureye.lib.shapefun;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import edu.uta.futureye.core.CoordinateTransform;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.function.AbstractMathFunc;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.FAxpb;
import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.basic.FX;
import edu.uta.futureye.function.intf.MathFunc;
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
public class SFQuadraticLocal2D extends AbstractMathFunc implements ScalarShapeFunction {
	private int funIndex;
	private MathFunc funCompose = null;
	private MathFunc funOuter = null;
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
		Map<String, MathFunc> fInners = new HashMap<String, MathFunc>(4);
		
		for(final String varName : varNames) {
			fInners.put(varName, new AbstractMathFunc(innerVarNames.toList()) {
				
				protected CoordinateTransform trans = new CoordinateTransform(2);
				
				public MathFunc diff(String var) {
					//Coordinate transform and Jacbian on element e
					List<MathFunc> funs = trans.getTransformFunction(
							//两种变换都可以，二次的要慢一些
							//trans.getTransformShapeFunctionByElement(e)
							trans.getTransformLinear2DShapeFunction(e)
								);
					trans.setTransformFunction(funs);
					
					MathFunc fx = funs.get(0);
					MathFunc fy = funs.get(1);
					
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

					MathFunc x_r = fx.diff("r");
					MathFunc x_s = fx.diff("s");
					MathFunc y_r = fy.diff("r");
					MathFunc y_s = fy.diff("s");
					
					//Function jac = (Function) trans.getJacobian2D();
					MathFunc jac = (MathFunc) e.getJacobin();
					
					if(varName.equals("r")) {
						if(var.equals("x"))
							return y_s.D(jac);
						if(var.equals("y"))
							return FC.C0.S(x_s.D(jac));
					} else if(varName.equals("s")) {
						if(var.equals("x"))
							return FC.C0.S(y_r.D(jac));
						if(var.equals("y"))
							return x_r.D(jac);
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

				@Override
				public double apply(Variable v) {
					// TODO Auto-generated method stub
					return 0;
				}
			});
		}
		
		MathFunc fr = new FX("r");
		MathFunc fs = new FX("s");
		MathFunc ft = //new FX("t");
				FC.C1.S(fr.A(fs));
				
				
		MathFunc f2rm1 = new FAxpb("r",2.0,-1.0);
		MathFunc f2sm1 = new FAxpb("s",2.0,-1.0);
		MathFunc f2tm1 = //new FAxpb("t",2.0,-1.0);
				FC.C1.S(FMath.linearCombination(2.0, fr, 2.0, fs));


//		Map<String, Function> fInner_t = new HashMap<String, Function>(4);
//		fInner_t.put("t", FOBasicDerivable.Minus(
//				new FConstant(1.0), 
//				FOBasicDerivable.Plus(fr, fs)
//				));
//		Function ftCompose = FOBasicDerivable.Compose(ft, fInner_t);
//		Function f2tm1Compose = FOBasicDerivable.Compose(f2tm1, fInner_t);
		MathFunc ftCompose = ft;
		MathFunc f2tm1Compose = f2tm1;
		
		if(funIndex == 0)
			funOuter = f2rm1.M(fr);
		else if(funIndex == 1)
			funOuter = f2sm1.M(fs);
		else if(funIndex == 2)
			funOuter = f2tm1Compose.M(ftCompose);
		else if(funIndex == 3)
			funOuter = FC.c(4.0).M(fr.M(fs));
		else if(funIndex == 4)
			funOuter = FC.c(4.0).M(fs.M(ftCompose));
		else if(funIndex == 5)
			funOuter = FC.c(4.0).M(fr.M(ftCompose));
		
		//设置funOuter的独立自变量名称
//		List<String> independentVarNames = new LinkedList<String>();
//		independentVarNames.add("r");
//		independentVarNames.add("s");
//		funOuter.setVarNames(independentVarNames);
		funOuter.setVarNames(varNames);
		//使用复合函数构造形函数
		funCompose = funOuter.compose(fInners);
		
	}

	@Override
	public void assignElement(Element e) {
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
	public MathFunc diff(String varName) {
		return funCompose.diff(varName);
	}

	@Override
	public double apply(Variable v) {
		return funCompose.apply(v);
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
			System.out.println(s.diff("r"));
		}
	}

	@Override
	public ObjList<String> innerVarNames() {
		return innerVarNames;
	}

}
