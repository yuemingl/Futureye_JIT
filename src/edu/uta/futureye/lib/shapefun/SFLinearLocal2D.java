package edu.uta.futureye.lib.shapefun;

import java.util.HashMap;
import java.util.Map;

import edu.uta.futureye.core.Element;
import edu.uta.futureye.function.AbstractFunction;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.function.intf.ScalarShapeFunction;
import edu.uta.futureye.util.FutureyeException;
import edu.uta.futureye.util.container.ObjList;
import edu.uta.futureye.util.container.VertexList;

/**
 * 三角形局部坐标，线性型函数
 *   Ni = N(r,s,t) = N( r(x,y), s(x,y), t(x,y) ), i=1,2,3
 *     N1 = r
 *     N2 = s
 *     N3 = t
 * where 
 *   r + s + t = 1
 * @author liuyueming
 *
 */
public class SFLinearLocal2D  extends AbstractFunction 
							  implements ScalarShapeFunction {
	protected int funIndex;
	private Function funOuter = null;
	private Function funCompose = null;
	protected ObjList<String> innerVarNames = null;
	
	protected Element e = null;
	private double area = -1.0;
	private double[] a = new double[3];
	private double[] b = new double[3];
	private double[] c = new double[3];
	
	private double coef = 1.0;
	
	class SF123 extends AbstractFunction {
		public SF123() {
			super(SFLinearLocal2D.this.varNames);
		}
		@Override
		public Function _d(String var) {
			if(varNames.get(funIndex).equals(var)) { 
				//d(N1)/dr = 1.0
				//d(N2)/ds = 1.0
				//d(N3)/dt = 1.0
				return new FC(1.0);
			} else if(funIndex == 2){ 
				//N3 = r = 1 - s - t, not free variable
				//d(N3)/ds = -1.0
				//d(N3)/dt = -1.0
				return new FC(-1.0);
			} else {
				return new FC(0.0);
			}
		}
		@Override
		public double value(Variable v) {
			return v.get(varNames.get(funIndex));
		}
		public String toString() {
			return varNames.get(funIndex);
		}
 	}
	
	/**
	 * 构造下列形函数中的一个：
	 * N1 = L1 = r
	 * N2 = L2 = s
	 * N3 = L3 = t
	 * @param funID = 1,2,3
	 * 
	 */
	public void Create(int funID, double coef) {
		funIndex = funID - 1;
		if( funID<1 || 3<funID ) {
			throw new FutureyeException("ERROR: funID should be 1,2 or 3.");
		}
		
		varNames.add("r");
		varNames.add("s");
		varNames.add("t");
		innerVarNames = new ObjList<String>("x","y");
		
		//复合函数
		Map<String, Function> fInners = new HashMap<String, Function>();
		
		final String varName = varNames.get(funIndex);
		fInners.put(varName, new AbstractFunction(innerVarNames.toList()) {	
			public Function _d(String var) {
				if(area < 0.0) {
					throw new FutureyeException("Check nodes order: area < 0.0");
				} else {
					if(varName.equals("r")) {//r对应三角形高h的负倒数
						if(var.equals("x"))
							return new FC(b[0] / (2 * area));
						if(var.equals("y"))
							return new FC(c[0] / (2 * area));
					} else if(varName.equals("s")) {//s对应三角形高h的负倒数
						if(var.equals("x"))
							return new FC(b[1] / (2 * area));
						if(var.equals("y"))
							return new FC(c[1] / (2 * area));
					} else if(varName.equals("t")) {//t对应三角形高h的负倒数
						if(var.equals("x"))
							return new FC(b[2] / (2 * area));
						if(var.equals("y"))
							return new FC(c[2] / (2 * area));
					}
				}
				return null;
			}
		});
		
		//使用复合函数构造形函数
		funOuter = new SF123();
		this.coef = coef;
		funCompose = FC.c(this.coef).M(
				funOuter.compose(fInners)
				); 
	}
	
	public SFLinearLocal2D(int funID) {
		this.Create(funID, 1.0);
	}
	
	public SFLinearLocal2D(int funID, double coef) {
		this.Create(funID, coef);
	}
	
	@Override
	public Function _d(String varName) {
		return funCompose._d(varName);
	}

	@Override
	public double value(Variable v) {
		return funCompose.value(v);
	}

	@Override
	public void asignElement(Element e) {
		this.e = e;
		//由node改为vertex，因为Element.adjustVerticeToCounterClockwise()结点顺序只调整了vertex
		VertexList vList = e.vertices();
		double x1 = vList.at(1).coord(1) , y1 =  vList.at(1).coord(2) ;
		double x2 = vList.at(2).coord(1) , y2 =  vList.at(2).coord(2) ;
		double x3 = vList.at(3).coord(1) , y3 =  vList.at(3).coord(2) ;
		
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
	}

	public String toString() {
		String varName = varNames.get(funIndex);
		return "N"+(funIndex+1)+"( "+varName+"(x,y) )="+funOuter.toString();
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
