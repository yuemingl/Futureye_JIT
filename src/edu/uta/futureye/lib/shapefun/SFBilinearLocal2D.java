package edu.uta.futureye.lib.shapefun;

import java.util.HashMap;
import java.util.Map;

import edu.uta.futureye.core.CoordinateTransform;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.function.MultiVarFunc;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.VariableArray;
import edu.uta.futureye.function.basic.FAxpb;
import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.function.intf.ScalarShapeFunction;
import edu.uta.futureye.util.Constant;
import edu.uta.futureye.util.FutureyeException;
import edu.uta.futureye.util.Utils;
import edu.uta.futureye.util.container.ObjList;
import edu.uta.futureye.util.container.VertexList;
import static edu.uta.futureye.function.FMath.*;

public class SFBilinearLocal2D extends MultiVarFunc implements ScalarShapeFunction {
	private int funIndex;
	private MathFunc funCompose = null;
	private MathFunc funOuter = null;
	private ObjList<String> innerVarNames = null;
	private double coef = 1.0;

	Element e;
	CoordinateTransform trans = new CoordinateTransform(2);
	MathFunc jac = null;
	MathFunc x_r = null;
	MathFunc x_s = null;
	MathFunc y_r = null;
	MathFunc y_s = null;
	
	
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
		
		varNames = new String[]{"r", "s"};
		innerVarNames = new ObjList<String>("x","y");
		
		//复合函数
		Map<String, MathFunc> fInners = new HashMap<String, MathFunc>(3);
		
		//r = r(x,y)
		//s = s(x,y)
		for(final String varName : varNames) {
			fInners.put(varName, new MultiVarFunc(varName, innerVarNames.toList()) {
				
/**
How to get derivative r_x, r_y, s_x, s_y:

f(x,y) = g(r,s)
f_x = g_r*r_x + g_s*s_x  ---(1)
f_y = g_r*r_y + g_s*s_y  ---(2)

for (1), let f=x and f=y we get tow equations, solve them:
(x_r x_s)   (r_x)   (1)
(y_r y_s) * (s_x) = (0)

similarly, for (2):
(x_r x_s)   (r_y)   (0)
(y_r y_s) * (s_y) = (1)

        (x_r x_s)
Let J = (y_r y_s)

from the above four equations, we have:
 (r_x r_y)
 (s_x s_y) = inv(J)
 */				
				//Derivatives: r_x, r_y, s_x, s_y
				public MathFunc diff(String var) {
					if(varName.equals("r")) {
						if(var.equals("x")) //r_x
							return y_s.D(jac);
						else //r_y
							return C0.S(x_s.D(jac));
					} else if(varName.equals("s")) {
						if(var.equals("x")) //s_x
							return C0.S(y_r.D(jac));
						else //s_y
							return x_r.D(jac);
					}
					return null;
				}
				
				@Override
				public String getExpr() {
					if(varName.equals("r"))
						return "r(x,y)";
					else
						return "s(x,y)";
				}
				
				@Override
				public String toString() {
					return getExpr();
				}
				
				@Override
				public double apply(double... args) {
					throw new RuntimeException("Unimplemented method");
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
		funCompose = FC.c(this.coef).M(funOuter.compose(fInners));
		/**
		 * The default active variable names of a composite function is the inner variable names.
		 * Shape function needs the outer variable names as the active variable names.
		 */
		funCompose.setActiveVarByNames(funOuter.getVarNames());
		funCompose.setArgIdx(Utils.getIndexMap(this.getVarNames()));
	}

	public SFBilinearLocal2D(int funID,double coef) {
		Create(funID,coef);
	}
	
	public SFBilinearLocal2D(int funID) {
		Create(funID,1.0);
	}

	public MathFunc diff(String varName) {
		//????this.funCompose.setActiveVarNames(this.getVarNames());
		return funCompose.diff(varName);
	}

	public double apply(Variable v) {
		return funCompose.apply(v);
	}
	
	@Override
	public double[] applyAll(VariableArray v, Map<Object,Object> cache) {
		return funCompose.applyAll(v,cache);
	}

	@Override
	public void assignElement(Element e) {
		this.e = e;
		VertexList vList = e.vertices();
	
//		//从_d()移动到这里，速度又快一倍
//		//Coordinate transform and Jacbian on element e
//		List<Function> funs = trans.getTransformFunction(
//				trans.getTransformLinear2DShapeFunction(e)
//				);
//		trans.setTransformFunction(funs);
//		Function fx = funs.get(0); //x=x(r,s)
//		Function fy = funs.get(1); //y=y(r,s)
//		x_r = fx._d("r");
//		x_s = fx._d("s");
//		y_r = fy._d("r");
//		y_s = fy._d("s");
		CoordinateTransform trans = e.getCoordTrans();
		if(trans == null) 
			throw new FutureyeException("call Element.updateJacobin() before calling assignElement()");
		MathFunc [] funs = trans.getJacobianMatrix();
		x_r = funs[0];
		x_s = funs[1];
		y_r = funs[2];
		y_s = funs[3];
		//用面积计算Jacobin，速度要快一倍
		double area = Utils.getRectangleArea(vList)/4.0;
		if(Math.abs(area)<Constant.eps) throw new FutureyeException();
		jac = FC.c(area);
	}

	public String getExpr() {
		if(this.coef != 1.0)
			return "N"+(funIndex+1)+"(r,s)";
		else
			return "N"+(funIndex+1)+"(r,s)";
	}
	
	public String toString() {
		if(this.coef != 1.0)
			return "N"+(funIndex+1)+"(r,s) = "+this.coef+"*"+funOuter.getExpr();
		else
			return "N"+(funIndex+1)+"(r,s) = "+funOuter.getExpr();
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

	@Override
	public double apply(double... args) {
		//no need. see MathFuncBase.compose()
		//this.funCompose.setActiveVarByNames(this.getVarNames());
		return this.funCompose.apply(args);
	}
}
