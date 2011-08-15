package edu.uta.futureye.lib.shapefun;

import java.util.HashMap;
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
import edu.uta.futureye.util.Utils;
import edu.uta.futureye.util.container.ObjList;
import edu.uta.futureye.util.container.VertexList;

public class SFBilinearLocal2D extends AbstractFunction implements ScalarShapeFunction {
	private int funIndex;
	private Function funCompose = null;
	private Function funOuter = null;
	private ObjList<String> innerVarNames = null;
	private double coef = 1.0;

	private Element e;
	private double jacFast = 0.0;
	
	
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
				
				public Function _d(String var) {
					//Coordinate transform and Jacbian on element e
					List<Function> funs = trans.getTransformFunction(
							trans.getTransformLinear2DShapeFunction(e)
							);
					trans.setTransformFunction(funs);
					
					Function fx = funs.get(0);
					Function fy = funs.get(1);
					
					Function x_r = fx._d("r");
					Function x_s = fx._d("s");
					Function y_r = fy._d("r");
					Function y_s = fy._d("s");
					
					//Function jac = trans.getJacobian2D();
					//Faster 快一倍
					Function jac = FC.c(jacFast);
					
					if(varName.equals("r")) {
						if(var.equals("x"))
							return y_s.D(jac);
						if(var.equals("y"))
							return FC.c0.S(x_s.D(jac));
					} else if(varName.equals("s")) {
						if(var.equals("x"))
							return FC.c0.S(y_r.D(jac));
						if(var.equals("y"))
							return x_r.D(jac);
					}
					return null;
				}
			});
		}
		
		if(funIndex == 0)
			funOuter = new FAxpb("r",-0.5,0.5).M(
					new FAxpb("s",-0.5,0.5));
		else if(funIndex == 1)
			funOuter = new FAxpb("r",0.5,0.5).M( 
					new FAxpb("s",-0.5,0.5));
		else if(funIndex == 2)
			funOuter = new FAxpb("r",0.5,0.5).M( 
					new FAxpb("s",0.5,0.5));
		else if(funIndex == 3)
			funOuter = new FAxpb("r",-0.5,0.5).M( 
					new FAxpb("s",0.5,0.5));
		
		//使用复合函数构造形函数
		this.coef = coef;
		funCompose = FC.c(this.coef).M(
					funOuter.compose(fInners)
				);
	}

	public SFBilinearLocal2D(int funID,double coef) {
		Create(funID,coef);
	}
	
	public SFBilinearLocal2D(int funID) {
		Create(funID,1.0);
	}

	public Function _d(String varName) {
		return funCompose._d(varName);
	}

	public double value(Variable v) {
		return funCompose.value(v);
	}

	@Override
	public void asignElement(Element e) {
		this.e = e;
		VertexList vList = e.vertices();
		//用面积计算Jacobin
		jacFast = Utils.getRectangleArea(vList)/4.0;
		
	}

	public String toString() {
		if(this.coef < 1.0)
			return "N"+(funIndex+1)+": "+this.coef+"*"+funOuter.toString();
		else
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
