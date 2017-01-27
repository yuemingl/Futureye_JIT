package edu.uta.futureye.lib.shapefun;

import java.util.HashMap;
import java.util.Map;

import edu.uta.futureye.core.Element;
import edu.uta.futureye.function.FMath;
import edu.uta.futureye.function.MultiVarFunc;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.FAxpb;
import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.basic.FX;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.function.intf.ScalarShapeFunction;
import edu.uta.futureye.util.container.ObjList;

/**
 * 2D Serendipity 局部坐标形函数
 * 
 * 
 * @author liuyueming
 *
 */
public class SFSerendipity2D extends MultiVarFunc implements ScalarShapeFunction {
	private int funIndex;
	private MathFunc funCompose = null;
	private MathFunc funOuter = null;
	private ObjList<String> innerVarNames = null;
	
	private MathFunc jac = null;
	private MathFunc x_r = null;
	private MathFunc x_s = null;
	private MathFunc y_r = null;
	private MathFunc y_s = null;
	
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
		
		varNames[0] = "r";
		varNames[1] = "s";
		innerVarNames = new ObjList<String>("x","y");
		
		//复合函数
		Map<String, MathFunc> fInners = new HashMap<String, MathFunc>(4);
		
		for(final String varName : varNames) {
			fInners.put(varName, new MultiVarFunc(varName, innerVarNames.toList()) {
				public MathFunc diff(String var) {
					if(varName.equals("r")) {
						if(var.equals("x"))
							return y_s.D(jac);
						if(var.equals("y"))
							return FMath.C0.S(x_s.D(jac));
					} else if(varName.equals("s")) {
						if(var.equals("x"))
							return FMath.C0.S(y_r.D(jac));
						if(var.equals("y"))
							return x_r.D(jac);
					}
					return null;
				}

				@Override
				public double apply(Variable v) {
					// TODO Auto-generated method stub
					return 0;
				}

				@Override
				public double apply(double... args) {
					// TODO Auto-generated method stub
					return 0;
				}
			});
		}
	
		MathFunc fx = new FX("r");
		MathFunc fy = new FX("s");
	
		MathFunc f1mx = new FAxpb("r",-1.0,1.0);
		MathFunc f1px = new FAxpb("r",1.0,1.0);
		MathFunc f1my = new FAxpb("s",-1.0,1.0);
		MathFunc f1py = new FAxpb("s",1.0,1.0);
		
		if(funIndex == 0)
			funOuter = FC.c(-0.25).M(f1mx).M(f1my).M(f1px.A(fy));
		else if(funIndex == 1)
			funOuter = FC.c(-0.25).M(f1px).M(f1my).M(f1mx.A(fy));
		else if(funIndex == 2)
			funOuter = FC.c(-0.25).M(f1px).M(f1py).M(f1mx.S(fy));
		else if(funIndex == 3)
			funOuter = FC.c(-0.25).M(f1mx).M(f1py).M(f1px.S(fy));
		else if(funIndex == 4)
			funOuter = FC.c(0.5).M(f1my).M(FMath.C1.S(fx.M(fx)));
		else if(funIndex == 5)
			funOuter = FC.c(0.5).M(f1px).M(FMath.C1.S(fy.M(fy)));
		else if(funIndex == 6)
			funOuter = FC.c(0.5).M(f1py).M(FMath.C1.S(fx.M(fx)));
		else if(funIndex == 7)
			funOuter = FC.c(0.5).M(f1mx).M(FMath.C1.S(fy.M(fy)));

		//使用复合函数构造形函数
		funCompose = funOuter.compose(fInners);		
	}
	
	@Override
	public void assignElement(Element e) {
		MathFunc[] funs = e.getCoordTrans().getJacobianMatrix();
		
		x_r = funs[0];
		x_s = funs[1];
		y_r = funs[2];
		y_s = funs[3];
		jac = e.getCoordTrans().getJacobian();
		
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
		return "N"+(funIndex+1)+"( r,s )="+funOuter.toString();
	}
	
	public static void main(String[] args){
		for(int i=1;i<=8;i++) {
			SFSerendipity2D s = new SFSerendipity2D(i);
			System.out.println(s);
			System.out.println(s.diff("r"));
			System.out.println(s.diff("s"));
		}
	}

	@Override
	public ObjList<String> innerVarNames() {
		return innerVarNames;
	}

	@Override
	public double apply(double... args) {
		// TODO Auto-generated method stub
		return 0;
	}
}
