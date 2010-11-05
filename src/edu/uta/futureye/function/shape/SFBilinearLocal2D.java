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
import edu.uta.futureye.function.intf.FunctionDerivable;
import edu.uta.futureye.function.intf.ShapeFunction;
import edu.uta.futureye.function.operator.FOBasicDerivable;

public class SFBilinearLocal2D implements ShapeFunction {
	private int funIndex;
	private FunctionDerivable funCompose = null;
	private FunctionDerivable funOuter = null;
	private List<String> varNames = new LinkedList<String>();
	private double coef = 1.0;

	private Element e;
	
	
	/**
	 * 构造下列形函数中的一个：
	 *  s
	 *  ^
	 *  |
	 *  |
	 * 
	 *  4-----3
	 *  |     |
	 *  |     |
	 *  1-----2  --> r
	 * -1  0  1
	 *
	 * N1 = (1-r)*(1-s)/4
	 * N2 = (1+r)*(1-s)/4
	 * N3 = (1+r)*(1+s)/4
	 * N4 = (1-r)*(1+s)/4
	 * @param funID = 1,2,3,4
	 * 
	 */	
	public void Create(int funID,double coef) {
		funIndex = funID - 1;
		if(funID<1 || funID>4) {
			System.out.println("ERROR: funID should be 1,2,3 or 4.");
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
							trans.getTransformLinear2DShapeFunction(e)
							);
					trans.setTransformFunction(funs);
					
					DerivativeIndicator di_r = new DerivativeIndicator("r1");
					DerivativeIndicator di_s = new DerivativeIndicator("s1");
					
					FunctionDerivable fx = funs.get(0);
					FunctionDerivable fy = funs.get(1);
					
					FunctionDerivable x_r = fx.derivative(di_r);
					FunctionDerivable x_s = fx.derivative(di_s);
					FunctionDerivable y_r = fy.derivative(di_r);
					FunctionDerivable y_s = fy.derivative(di_s);
					FunctionDerivable jac = trans.getJacobian2D();
					
					//TODO ??? varName.equals("r")
					if(varName == "r") {
						if(di.get().get("x") != null)
							return FOBasicDerivable.Divi(y_s, jac);
						if(di.get().get("y") != null)
							return FOBasicDerivable.Minus(new FConstant(0.0),FOBasicDerivable.Divi(x_s, jac));
					} else if(varName == "s") {
						if(di.get().get("x") != null)
							return FOBasicDerivable.Minus(new FConstant(0.0),FOBasicDerivable.Divi(y_r, jac));
						if(di.get().get("y") != null)
							return FOBasicDerivable.Divi(x_r, jac);
					}
					return null;
				}
			});
		}
		
		if(funIndex == 0)
			funOuter = FOBasicDerivable.Mult(new FAxpb("r",-0.5,0.5), 
					new FAxpb("s",-0.5,0.5));
		else if(funIndex == 1)
			funOuter = FOBasicDerivable.Mult(new FAxpb("r",0.5,0.5), 
					new FAxpb("s",-0.5,0.5));
		else if(funIndex == 2)
			funOuter = FOBasicDerivable.Mult(new FAxpb("r",0.5,0.5), 
					new FAxpb("s",0.5,0.5));
		else if(funIndex == 3)
			funOuter = FOBasicDerivable.Mult(new FAxpb("r",-0.5,0.5), 
					new FAxpb("s",0.5,0.5));
		
		//使用复合函数构造形函数
		this.coef = coef;
		funCompose = FOBasicDerivable.Mult(new FConstant(this.coef), 
				FOBasicDerivable.Compose(funOuter, fInners));
	}

	public SFBilinearLocal2D(int funID,double coef) {
		Create(funID,coef);
	}
	
	public SFBilinearLocal2D(int funID) {
		Create(funID,1.0);
	}

	public FunctionDerivable derivative(DerivativeIndicator di) {
		return funCompose.derivative(di);
	}

	public double value(Variable v) {
		return funCompose.value(v);
	}

	@Override
	public void asignElement(Element e) {
		this.e = e;
	}

	@Override
	public void setVarNames(List<String> varNames) {
		this.varNames = varNames;
	}

	@Override
	public List<String> varNames() {
		return varNames;
	}
	
	public String toString() {
		return "N"+(funIndex+1)+": "+funOuter.toString();
	}

	ShapeFunction sf1d1 = new SFLinearLocal1D(1);
	ShapeFunction sf1d2 = new SFLinearLocal1D(2);
	@Override
	public ShapeFunction restrictTo(int funIndex) {
		if(funIndex == 1) return sf1d1;
		else return sf1d2;
	}
}
