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
import edu.uta.futureye.function.intf.FunctionDerivable;
import edu.uta.futureye.function.intf.ShapeFunction;
import edu.uta.futureye.function.operator.FOBasicDerivable;

public class SFSerendipity2D implements ShapeFunction {
	private int funIndex;
	private FunctionDerivable funCompose = null;
	private FunctionDerivable funOuter = null;
	private List<String> varNames = new LinkedList<String>();
	private Element e = null;
	/**
	 * 构造下列形函数中的一个：
	 *  s
	 *  ^
	 *  |
	 *  |
	 * 
	 *  4--7--3
	 *  |     |
	 *  8     6
	 *  |     |
	 *  1--5--2  --> r
	 * -1  0  1
	 *
	 *for i=1,2,3,4
	 * Ni = (1+r0)*(1+s0)*(r0+s0-1)/4
	 * 
	 *for i=5,6,7,8
	 * Ni = (1-r^2)*(1+s0), when ri=0
	 * Ni = (1+r0)*(1-s^2), when si=0
	 * 
	 *where
	 * r0=r*ri
	 * s0=s*si
	 * 
	 * @param funID = 1,...,8
	 * 
	 */	
	public SFSerendipity2D(int funID) {
		funIndex = funID - 1;
		if(funID<1 || funID>8) {
			System.out.println("ERROR: funID should be 1,...,8.");
			return;
		}
		
		varNames.add("r");
		varNames.add("s");
		
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
							trans.getTransformShapeFunctionByElement(e)
								);
					trans.setTransformFunction(funs);
					
					DerivativeIndicator di_r = new DerivativeIndicator("r1");
					DerivativeIndicator di_s = new DerivativeIndicator("s1");
					
					FunctionDerivable x_r = funs.get(0).derivative(di_r);
					FunctionDerivable x_s = funs.get(0).derivative(di_s);
					FunctionDerivable y_r = funs.get(1).derivative(di_r);
					FunctionDerivable y_s = funs.get(1).derivative(di_s);
					FunctionDerivable jac = trans.getJacobian2D();
					
					//TODO !!! 数值积分点的选取问题
//					Variable v = new Variable();
//					v.set("r", 0.0);
//					v.set("s", 0.0);
//					System.out.println(x_r.value(v));
//					System.out.println(x_s.value(v));
//					System.out.println(y_r.value(v));
//					System.out.println(y_s.value(v));
//					System.out.println(jac.value(v));
//					0.0
//					0.33333333349999994
//					-0.33333333349999994
//					0.0
//					0.11111111122222218
					
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
					
					return null;
				}
			});
		}
	
		FunctionDerivable fx = new FX("r");
		FunctionDerivable fy = new FX("s");
	
		FunctionDerivable f1mx = new FAxpb("r",-1.0,1.0);
		FunctionDerivable f1px = new FAxpb("r",1.0,1.0);
		FunctionDerivable f1my = new FAxpb("s",-1.0,1.0);
		FunctionDerivable f1py = new FAxpb("s",1.0,1.0);
		
		if(funIndex == 0)
			funOuter = FOBasicDerivable.Mult(
					FOBasicDerivable.Mult(
					FOBasicDerivable.Mult(new FConstant(-0.25), f1mx),f1my),
					FOBasicDerivable.Plus(f1px, fy));
		else if(funIndex == 1)
			funOuter = FOBasicDerivable.Mult(
					FOBasicDerivable.Mult(
					FOBasicDerivable.Mult(new FConstant(-0.25), f1px),f1my),
					FOBasicDerivable.Plus(f1mx, fy));
		else if(funIndex == 2)
			funOuter = FOBasicDerivable.Mult(
					FOBasicDerivable.Mult(
					FOBasicDerivable.Mult(new FConstant(-0.25), f1px),f1py),
					FOBasicDerivable.Minus(f1mx, fy));
		else if(funIndex == 3)
			funOuter = FOBasicDerivable.Mult(
					FOBasicDerivable.Mult(
					FOBasicDerivable.Mult(new FConstant(-0.25), f1mx),f1py),
					FOBasicDerivable.Minus(f1px, fy));
		else if(funIndex == 4)
			funOuter = FOBasicDerivable.Mult(
					FOBasicDerivable.Mult(new FConstant(0.5), f1my),
					FOBasicDerivable.Minus(new FConstant(1.0), 
							FOBasicDerivable.Mult(fx,fx)));
		else if(funIndex == 5)
			funOuter = FOBasicDerivable.Mult(
					FOBasicDerivable.Mult(new FConstant(0.5), f1px),
					FOBasicDerivable.Minus(new FConstant(1.0), 
							FOBasicDerivable.Mult(fy,fy)));
		else if(funIndex == 6)
			funOuter = FOBasicDerivable.Mult(
					FOBasicDerivable.Mult(new FConstant(0.5), f1py),
					FOBasicDerivable.Minus(new FConstant(1.0), 
							FOBasicDerivable.Mult(fx,fx)));
		else if(funIndex == 7)
			funOuter = FOBasicDerivable.Mult(
					FOBasicDerivable.Mult(new FConstant(0.5), f1mx),
					FOBasicDerivable.Minus(new FConstant(1.0), 
							FOBasicDerivable.Mult(fy,fy)));

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
		return "N"+(funIndex+1)+"( r,s )="+funOuter.toString();
	}
	
	public static void main(String[] args){
		for(int i=1;i<=8;i++) {
			SFSerendipity2D s = new SFSerendipity2D(i);
			System.out.println(s);
			DerivativeIndicator di = new DerivativeIndicator();
			di.set("r", 1);
			System.out.println(s.derivative(di));
			di = new DerivativeIndicator();
			di.set("s", 1);
			System.out.println(s.derivative(di));
		}
	}
}
