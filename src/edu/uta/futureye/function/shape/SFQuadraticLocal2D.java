package edu.uta.futureye.function.shape;

import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import edu.uta.futureye.core.CoordinateTransform;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.function.DerivativeIndicator;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.FAbstract;
import edu.uta.futureye.function.basic.FAxpb;
import edu.uta.futureye.function.basic.FConstant;
import edu.uta.futureye.function.basic.FX;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.function.intf.FunctionDerivable;
import edu.uta.futureye.function.intf.ShapeFunction;
import edu.uta.futureye.function.operator.FOBasicDerivable;

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
public class SFQuadraticLocal2D implements ShapeFunction {
	private int funIndex;
	private FunctionDerivable funCompose = null;
	private FunctionDerivable funOuter = null;
	private List<String> varNames = new LinkedList<String>();

	private Element e = null;

	public SFQuadraticLocal2D(int funID) {
		funIndex = funID - 1;
		if(funID<1 || funID>6) {
			System.out.println("ERROR: funID should be 1~6.");
			return;
		}
		
		varNames.add("r");
		varNames.add("s");
		//varNames.add("t");
		
		//复合函数
		Map<String, FunctionDerivable> fInners = new HashMap<String, FunctionDerivable>(4);
		List<String> varNamesInner = new LinkedList<String>();
		varNamesInner.add("x");
		varNamesInner.add("y");
		
		for(final String varName : varNames) {
			fInners.put(varName, new FAbstract(varNamesInner) {
				
				protected CoordinateTransform trans = new CoordinateTransform(2);
				
				public FunctionDerivable derivative(DerivativeIndicator di) {
					//Coordinate transform and Jacbian on element e
					List<FunctionDerivable> funs = trans.getTransformFunction(
							//两种变换都可以，二次的要慢一些
							//trans.getTransformShapeFunctionByElement(e)
							trans.getTransformLinear2DShapeFunction(e)
								);
					trans.setTransformFunction(funs);
					
					DerivativeIndicator di_r = new DerivativeIndicator("r1");
					DerivativeIndicator di_s = new DerivativeIndicator("s1");
					
					FunctionDerivable fx = funs.get(0);
					FunctionDerivable fy = funs.get(1);
					
//					FunctionDerivable fr = new FX("r");
//					FunctionDerivable fs = new FX("s");
//					Map<String, FunctionDerivable> fInner_t = new HashMap<String, FunctionDerivable>(4);
//					fInner_t.put("t", FOBasicDerivable.Minus(
//							new FConstant(1.0), 
//							FOBasicDerivable.Plus(fr, fs)
//							));
//					FunctionDerivable fxCompose = FOBasicDerivable.Compose(fx, fInner_t);
//					FunctionDerivable fyCompose = FOBasicDerivable.Compose(fy, fInner_t);
//					
//					List<String> independentVarNames = new LinkedList<String>();
//					independentVarNames.add("r");
//					independentVarNames.add("s");
//					fxCompose.setVarNames(independentVarNames);
//					fyCompose.setVarNames(independentVarNames);

					FunctionDerivable x_r = fx.derivative(di_r);
					FunctionDerivable x_s = fx.derivative(di_s);
					FunctionDerivable y_r = fy.derivative(di_r);
					FunctionDerivable y_s = fy.derivative(di_s);
					
					//FunctionDerivable jac = (Function) trans.getJacobian2D();
					FunctionDerivable jac = (FunctionDerivable) e.getJacobin();
					
					if(varName.equals("r")) {
						if(di.get().get("x") != null)
							return FOBasicDerivable.Divi(y_s, jac);
						if(di.get().get("y") != null)
							return FOBasicDerivable.Minus(new FConstant(0.0),FOBasicDerivable.Divi(x_s, jac));
					} else if(varName.equals("s")) {
						if(di.get().get("x") != null)
							return FOBasicDerivable.Minus(new FConstant(0.0),FOBasicDerivable.Divi(y_r, jac));
						if(di.get().get("y") != null)
							return FOBasicDerivable.Divi(x_r, jac);
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
		
		FunctionDerivable fr = new FX("r");
		FunctionDerivable fs = new FX("s");
		FunctionDerivable ft = //new FX("t");
		FOBasicDerivable.Minus(
				new FConstant(1.0), 
				FOBasicDerivable.Plus(fr, fs));
				
				
		FunctionDerivable f2rm1 = new FAxpb("r",2.0,-1.0);
		FunctionDerivable f2sm1 = new FAxpb("s",2.0,-1.0);
		FunctionDerivable f2tm1 = //new FAxpb("t",2.0,-1.0);
			FOBasicDerivable.Minus(
					new FConstant(1.0), 
					FOBasicDerivable.LinearCombination(2.0, fr, 2.0, fs));


//		Map<String, FunctionDerivable> fInner_t = new HashMap<String, FunctionDerivable>(4);
//		fInner_t.put("t", FOBasicDerivable.Minus(
//				new FConstant(1.0), 
//				FOBasicDerivable.Plus(fr, fs)
//				));
//		FunctionDerivable ftCompose = FOBasicDerivable.Compose(ft, fInner_t);
//		FunctionDerivable f2tm1Compose = FOBasicDerivable.Compose(f2tm1, fInner_t);
		FunctionDerivable ftCompose = ft;
		FunctionDerivable f2tm1Compose = f2tm1;
		
		if(funIndex == 0)
			funOuter = FOBasicDerivable.Mult(f2rm1, fr);
		else if(funIndex == 1)
			funOuter = FOBasicDerivable.Mult(f2sm1, fs);
		else if(funIndex == 2)
			funOuter = FOBasicDerivable.Mult(f2tm1Compose, ftCompose);
		else if(funIndex == 3)
			funOuter = FOBasicDerivable.Mult(new FConstant(4.0), FOBasicDerivable.Mult(fr, fs));
		else if(funIndex == 4)
			funOuter = FOBasicDerivable.Mult(new FConstant(4.0), FOBasicDerivable.Mult(fs, ftCompose));
		else if(funIndex == 5)
			funOuter = FOBasicDerivable.Mult(new FConstant(4.0), FOBasicDerivable.Mult(fr, ftCompose));
		
		//设置funOuter的独立自变量名称
//		List<String> independentVarNames = new LinkedList<String>();
//		independentVarNames.add("r");
//		independentVarNames.add("s");
//		funOuter.setVarNames(independentVarNames);
		funOuter.setVarNames(varNames);
		//使用复合函数构造形函数
		funCompose = FOBasicDerivable.Compose(funOuter, fInners);
		
	}

	@Override
	public void asignElement(Element e) {
		this.e = e;
	}

	//TODO ??? 应该采用一维二次型函数
	ShapeFunction sf1d1 = new SFLinearLocal1D(1);
	ShapeFunction sf1d2 = new SFLinearLocal1D(2);
	@Override
	public ShapeFunction restrictTo(int funIndex) {
		if(funIndex == 1) return sf1d1;
		else return sf1d2;
	}

	@Override
	public FunctionDerivable derivative(DerivativeIndicator di) {
		return funCompose.derivative(di);
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
			DerivativeIndicator di = new DerivativeIndicator();
			di.set("r", 1);
			System.out.println(s.derivative(di));			
		}
	}

}
