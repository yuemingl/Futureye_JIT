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
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.function.intf.ScalarShapeFunction;
import edu.uta.futureye.function.operator.FMath;
import edu.uta.futureye.util.list.ObjList;

public class SFBilinearLocal2D extends AbstractFunction implements ScalarShapeFunction {
	private int funIndex;
	private Function funCompose = null;
	private Function funOuter = null;
	private List<String> varNames = new LinkedList<String>();
	private ObjList<String> innerVarNames = null;
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
		innerVarNames = new ObjList<String>("x","y");
		
		//复合函数
		Map<String, Function> fInners = new HashMap<String, Function>(4);
		
		for(final String varName : varNames) {
			fInners.put(varName, new AbstractFunction(innerVarNames.toList()) {
				
				protected CoordinateTransform trans = new CoordinateTransform(2);
				
				public Function d(String var) {
					//Coordinate transform and Jacbian on element e
					List<Function> funs = trans.getTransformFunction(
							trans.getTransformLinear2DShapeFunction(e)
							);
					trans.setTransformFunction(funs);
					
					Function fx = funs.get(0);
					Function fy = funs.get(1);
					
					Function x_r = fx.d("r");
					Function x_s = fx.d("s");
					Function y_r = fy.d("r");
					Function y_s = fy.d("s");
					Function jac = trans.getJacobian2D();
					
					if(varName.equals("r")) {
						if(var.equals("x"))
							return FMath.Divi(y_s, jac);
						if(var.equals("y"))
							//return FOBasicDerivable.Divi(x_s, jac);
							return FMath.Minus(new FC(0.0),FMath.Divi(x_s, jac));
					} else if(varName.equals("s")) {
						if(var.equals("x"))
							//return FOBasicDerivable.Divi(y_r, jac);
							return FMath.Minus(new FC(0.0),FMath.Divi(y_r, jac));
						if(var.equals("y"))
							return FMath.Divi(x_r, jac);
					}
					return null;
				}
			});
		}
		
		if(funIndex == 0)
			funOuter = FMath.Mult(new FAxpb("r",-0.5,0.5), 
					new FAxpb("s",-0.5,0.5));
		else if(funIndex == 1)
			funOuter = FMath.Mult(new FAxpb("r",0.5,0.5), 
					new FAxpb("s",-0.5,0.5));
		else if(funIndex == 2)
			funOuter = FMath.Mult(new FAxpb("r",0.5,0.5), 
					new FAxpb("s",0.5,0.5));
		else if(funIndex == 3)
			funOuter = FMath.Mult(new FAxpb("r",-0.5,0.5), 
					new FAxpb("s",0.5,0.5));
		
		//使用复合函数构造形函数
		this.coef = coef;
		funCompose = FMath.Mult(new FC(this.coef), 
				FMath.Compose(funOuter, fInners));
	}

	public SFBilinearLocal2D(int funID,double coef) {
		Create(funID,coef);
	}
	
	public SFBilinearLocal2D(int funID) {
		Create(funID,1.0);
	}

	public Function d(String varName) {
		return funCompose.d(varName);
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

	ScalarShapeFunction sf1d1 = new SFLinearLocal1D(1);
	ScalarShapeFunction sf1d2 = new SFLinearLocal1D(2);
	@Override
	public ScalarShapeFunction restrictTo(int funIndex) {
		if(funIndex == 1) return sf1d1;
		else return sf1d2;
	}

	@Override
	public ObjList<String> innerVarNames() {
		return innerVarNames;
	}
}
