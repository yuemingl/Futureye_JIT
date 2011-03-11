package edu.uta.futureye.lib.shapefun;

import java.util.LinkedList;
import java.util.List;

import edu.uta.futureye.core.Element;
import edu.uta.futureye.function.AbstractFunction;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.FAxpb;
import edu.uta.futureye.function.basic.FX;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.function.intf.ScalarShapeFunction;
import edu.uta.futureye.function.operator.FMath;
import edu.uta.futureye.util.FutureyeException;
import edu.uta.futureye.util.list.ObjList;

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
public class SFQuadraticLocal2DFast extends AbstractFunction implements ScalarShapeFunction {
	private int funIndex;
	private Function funOuter = null;
	private ObjList<String> innerVarNames = null;

	private double area = -1.0;
	private double[] a = new double[3];
	private double[] b = new double[3];
	private double[] c = new double[3];
	private double[] dlx = new double[3];
	private double[] dly = new double[3];

	/**
	 * Real shape function N1-N6
	 * @author liuyueming
	 *
	 */
	protected class SF123456 extends AbstractFunction {
		int funIndex;
		public SF123456(int funIndex) {
			super(SFQuadraticLocal2DFast.this.varNames);
			this.funIndex = funIndex;
		}
		@Override
		public Function d(String var) {
			if(area < 0.0) {
				FutureyeException e = new FutureyeException("SFLinearLocal2D: area < 0.0");
				e.printStackTrace();
				return null;
			}
			Function fr = new FX("r");
			Function fs = new FX("s");
			Function ft = new FX("t");
			if(var.equals("x")) {
				if(funIndex == 0)
					return new FAxpb("r",dlx[0]*4.0,-dlx[0]);
				else if(funIndex == 1)
					return new FAxpb("s",dlx[1]*4.0,-dlx[1]);
				else if(funIndex == 2)
					return new FAxpb("t",dlx[2]*4.0,-dlx[2]);
				else if(funIndex == 3)
					return FMath.LinearCombination(4.0*dlx[0], fs, 4.0*dlx[1], fr);
				else if(funIndex == 4)
					return FMath.LinearCombination(4.0*dlx[1], ft, 4.0*dlx[2], fs);
				else if(funIndex == 5)
					return FMath.LinearCombination(4.0*dlx[0], ft, 4.0*dlx[2], fr);		
				else {
					FutureyeException e = new FutureyeException("Error: SF123456.derivative(x)");
					e.printStackTrace();
				}
			} else if(var.equals("y")) {
				if(funIndex == 0)
					return new FAxpb("r",dly[0]*4.0,-dly[0]);
				else if(funIndex == 1)
					return new FAxpb("s",dly[1]*4.0,-dly[1]);
				else if(funIndex == 2)
					return new FAxpb("t",dly[2]*4.0,-dly[2]);
				else if(funIndex == 3)
					return FMath.LinearCombination(4.0*dly[0], fs, 4.0*dly[1], fr);
				else if(funIndex == 4)
					return FMath.LinearCombination(4.0*dly[1], ft, 4.0*dly[2], fs);
				else if(funIndex == 5)
					return FMath.LinearCombination(4.0*dly[0], ft, 4.0*dly[2], fr);		
				else {
					FutureyeException e = new FutureyeException("Error: SF123456.derivative(y)");
					e.printStackTrace();
				}
			} else {
				FutureyeException e = new FutureyeException("Error: SF123456.derivative()");
				e.printStackTrace();				
			}
			return null;
		}
		
		/**
		 * N1 = (2*r-1)*r
		 * N2 = (2*s-1)*s
		 * N3 = (2*t-1)*t
		 * N4 = 4*r*s
		 * N5 = 4*s*t
		 * N6 = 4*r*t
		 */
		@Override
		public double value(Variable v) {
			double r = v.get("r");
			double s = v.get("s");
			double t = v.get("t");
			if(funIndex == 0)
				return (2.0*r-1.0)*r;
			else if(funIndex == 1)
				return (2.0*s-1.0)*s;
			else if(funIndex == 2)
				return (2.0*t-1.0)*t;
			else if(funIndex == 3)
				return 4.0*r*s;
			else if(funIndex == 4)
				return 4.0*s*t;
			else if(funIndex == 5)
				return 4.0*r*t;
			else {
				FutureyeException e = new FutureyeException("Error: SF123456.value()");
				e.printStackTrace();
			}
			return 0.0;
		}
		
		public String toString() {
			String[] s = {
					 "(2*r-1)*r",
					 "(2*s-1)*s",
					 "(2*t-1)*t",
					 "4*r*s",
					 "4*s*t",
					 "4*r*t"
			};
			return s[funIndex];
		}
 	}	
	
	
	public SFQuadraticLocal2DFast(int funID) {
		funIndex = funID - 1;
		if(funID<1 || funID>6) {
			System.out.println("ERROR: funID should be 1~6.");
			return;
		}
		
		varNames.add("r");
		varNames.add("s");
		varNames.add("t");
		innerVarNames = new ObjList<String>("x","y");
		
		funOuter = new SF123456(funIndex);
	}

	@Override
	public void asignElement(Element e) {
		
		double x1 = e.nodes.at(1).coord(1) , y1 =  e.nodes.at(1).coord(2) ;
		double x2 = e.nodes.at(2).coord(1) , y2 =  e.nodes.at(2).coord(2) ;
		double x3 = e.nodes.at(3).coord(1) , y3 =  e.nodes.at(3).coord(2) ;
		area = ( (x2*y3 - x3*y2) - (x1*y3 - x3*y1) + (x1*y2 - x2*y1) ) / 2.0;
		a[0] = x2*y3 - x3*y2;
		b[0] = y2 - y3;
		c[0] = x3 - x2;
		a[1] = x3*y1 - x1*y3;
		b[1] = y3 - y1;
		c[1] = x1 - x3;
		a[2] = x1*y2 - x2*y1;
		b[2] = y1 - y2;
		c[2] = x2 - x1;
		
		for(int i=0;i<3;i++) {
			dlx[i] = b[i]/(2.0*area);
			dly[i] = c[i]/(2.0*area);
		}
	}

	//TODO ??? 应该采用一维二次型函数
	ScalarShapeFunction sf1d1 = new SFQuadraticLocal1D(1);
	ScalarShapeFunction sf1d2 = new SFQuadraticLocal1D(2);
	ScalarShapeFunction sf1d3 = new SFQuadraticLocal1D(3);
	@Override
	public ScalarShapeFunction restrictTo(int funIndex) {
		if(funIndex == 1) return sf1d1;
		else if(funIndex == 2) return sf1d2;
		else if(funIndex == 3) return sf1d3;
		else {
			FutureyeException e = new FutureyeException("Error: SF123456.restrictTo()");
			e.printStackTrace();
			return null;
		}
	}

	@Override
	public Function d(String varName) {
		return funOuter.d(varName);
	}

	@Override
	public double value(Variable v) {
		return funOuter.value(v);
	}

	public String toString() {
		return "N"+(funIndex+1)+"( r,s,t )="+funOuter.toString();
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		for(int i=1;i<=6;i++) {
			SFQuadraticLocal2DFast sf = new SFQuadraticLocal2DFast(i);
			SFQuadraticLocal2D s = new SFQuadraticLocal2D(i);
			System.out.println(sf);
			System.out.println(s);
			Variable var = new Variable();
			var.set("r", 0.6); var.set("s", 0.2); var.set("t", 0.2);
			System.out.println(sf.value(var));
			System.out.println(s.value(var));
			//DerivativeIndicator di = new DerivativeIndicator();
			//di.set("x", 1);
			//System.out.println(s.derivative(di));			
		}
	}

	@Override
	public ObjList<String> innerVarNames() {
		return innerVarNames;
	}

}
