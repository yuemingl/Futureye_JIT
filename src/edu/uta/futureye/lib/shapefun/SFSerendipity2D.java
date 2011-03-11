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
import edu.uta.futureye.util.list.ObjList;

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
	private List<String> varNames = new LinkedList<String>();
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
				
				public Function d(String var) {
					//Coordinate transform and Jacbian on element e
					List<Function> funs = trans.getTransformFunction(
							trans.getTransformLinear2DShapeFunction(e)
							//TODO 下面的调用有问题
							//trans.getTransformShapeFunctionByElement(e)
								);
					trans.setTransformFunction(funs);
					
					Function x_r = funs.get(0).d("r");
					Function x_s = funs.get(0).d("s");
					Function y_r = funs.get(1).d("r");
					Function y_s = funs.get(1).d("s");
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
							return FMath.Divi(y_s, jac);
						if(var.equals("y"))
							return FMath.Minus(new FC(0.0),FMath.Divi(x_s, jac));
					} else if(varName.equals("s")) {
						if(var.equals("x"))
							return FMath.Minus(new FC(0.0),FMath.Divi(y_r, jac));
						if(var.equals("y"))
							return FMath.Divi(x_r, jac);
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
			funOuter = FMath.Mult(
					FMath.Mult(
					FMath.Mult(new FC(-0.25), f1mx),f1my),
					FMath.Plus(f1px, fy));
		else if(funIndex == 1)
			funOuter = FMath.Mult(
					FMath.Mult(
					FMath.Mult(new FC(-0.25), f1px),f1my),
					FMath.Plus(f1mx, fy));
		else if(funIndex == 2)
			funOuter = FMath.Mult(
					FMath.Mult(
					FMath.Mult(new FC(-0.25), f1px),f1py),
					FMath.Minus(f1mx, fy));
		else if(funIndex == 3)
			funOuter = FMath.Mult(
					FMath.Mult(
					FMath.Mult(new FC(-0.25), f1mx),f1py),
					FMath.Minus(f1px, fy));
		else if(funIndex == 4)
			funOuter = FMath.Mult(
					FMath.Mult(new FC(0.5), f1my),
					FMath.Minus(new FC(1.0), 
							FMath.Mult(fx,fx)));
		else if(funIndex == 5)
			funOuter = FMath.Mult(
					FMath.Mult(new FC(0.5), f1px),
					FMath.Minus(new FC(1.0), 
							FMath.Mult(fy,fy)));
		else if(funIndex == 6)
			funOuter = FMath.Mult(
					FMath.Mult(new FC(0.5), f1py),
					FMath.Minus(new FC(1.0), 
							FMath.Mult(fx,fx)));
		else if(funIndex == 7)
			funOuter = FMath.Mult(
					FMath.Mult(new FC(0.5), f1mx),
					FMath.Minus(new FC(1.0), 
							FMath.Mult(fy,fy)));

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
	public Function d(String varName) {
		return funCompose.d(varName);
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
			System.out.println(s.d("r"));
			System.out.println(s.d("s"));
		}
	}

	@Override
	public ObjList<String> innerVarNames() {
		return innerVarNames;
	}
}
