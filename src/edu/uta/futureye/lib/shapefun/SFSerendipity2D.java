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
import edu.uta.futureye.function.basic.FX;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.function.intf.ScalarShapeFunction;
import edu.uta.futureye.util.container.ObjList;

/**
 * 2D Serendipity 局部坐标形函数
 * 
 * 
 * @author liuyueming
 *
 */
public class SFSerendipity2D extends AbstractFunction implements ScalarShapeFunction {
	private int funIndex;
	private Function funCompose = null;
	private Function funOuter = null;
	private ObjList<String> innerVarNames = null;
	
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
							//TODO 下面的调用有问题
							//trans.getTransformShapeFunctionByElement(e)
								);
					trans.setTransformFunction(funs);
					
					Function x_r = funs.get(0)._d("r");
					Function x_s = funs.get(0)._d("s");
					Function y_r = funs.get(1)._d("r");
					Function y_s = funs.get(1)._d("s");
					Function jac = trans.getJacobian2D();
					
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
	
		Function fx = new FX("r");
		Function fy = new FX("s");
	
		Function f1mx = new FAxpb("r",-1.0,1.0);
		Function f1px = new FAxpb("r",1.0,1.0);
		Function f1my = new FAxpb("s",-1.0,1.0);
		Function f1py = new FAxpb("s",1.0,1.0);
		
		if(funIndex == 0)
			funOuter = FC.c(-0.25).M(f1mx).M(f1my).M(f1px.A(fy));
		else if(funIndex == 1)
			funOuter = FC.c(-0.25).M(f1px).M(f1my).M(f1mx.A(fy));
		else if(funIndex == 2)
			funOuter = FC.c(-0.25).M(f1px).M(f1py).M(f1mx.S(fy));
		else if(funIndex == 3)
			funOuter = FC.c(-0.25).M(f1mx).M(f1py).M(f1px.S(fy));
		else if(funIndex == 4)
			funOuter = FC.c(0.5).M(f1my).M(FC.c1.S(fx.M(fx)));
		else if(funIndex == 5)
			funOuter = FC.c(0.5).M(f1px).M(FC.c1.S(fy.M(fy)));
		else if(funIndex == 6)
			funOuter = FC.c(0.5).M(f1py).M(FC.c1.S(fx.M(fx)));
		else if(funIndex == 7)
			funOuter = FC.c(0.5).M(f1mx).M(FC.c1.S(fy.M(fy)));

		//使用复合函数构造形函数
		funCompose = funOuter.compose(fInners);		
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
	public double value(Variable v) {
		return funCompose.value(v);
	}

	public String toString() {
		return "N"+(funIndex+1)+"( r,s )="+funOuter.toString();
	}
	
	public static void main(String[] args){
		for(int i=1;i<=8;i++) {
			SFSerendipity2D s = new SFSerendipity2D(i);
			System.out.println(s);
			System.out.println(s._d("r"));
			System.out.println(s._d("s"));
		}
	}

	@Override
	public ObjList<String> innerVarNames() {
		return innerVarNames;
	}
}
