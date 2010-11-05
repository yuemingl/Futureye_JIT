package edu.uta.futureye.function.shape;

import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Vertex;
import edu.uta.futureye.function.DerivativeIndicator;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.FAbstract;
import edu.uta.futureye.function.basic.FConstant;
import edu.uta.futureye.function.intf.FunctionDerivable;
import edu.uta.futureye.function.intf.ShapeFunction;
import edu.uta.futureye.function.operator.FOBasicDerivable;
import edu.uta.futureye.util.FutureEyeException;

/**
 * 三角形局部坐标，线性型函数
 * N = N(r,s,t) = N( r(x,y), s(x,y), t(x,y) )
 * N1 = r
 * N2 = s
 * N3 = t
 * @author liuyueming
 *
 */
public class SFLinearLocal2D implements ShapeFunction {
	private int funIndex;
	private FunctionDerivable funCompose = null;
	private FunctionDerivable funOuter = null;
	private List<String> sfVarNames = new LinkedList<String>();
	
	protected Element e = null;
	private double area = -1.0;
	private double[] a = new double[3];
	private double[] b = new double[3];
	private double[] c = new double[3];
	
	private double coef = 1.0;
	
	protected class SF123 extends FAbstract {
		int funIndex;
		public SF123(int funIndex) {
			super(sfVarNames);
			this.funIndex = funIndex;
		}
		@Override
		public FunctionDerivable derivative(DerivativeIndicator di) {
			//2阶以上导数
			if( di.get().values().iterator().next() > 1)
				return new FConstant(0.0);
			else {
				String ss = sfVarNames.get(funIndex);
				Integer degree = di.get().get(ss);
				if(degree != null) {
					return new FConstant(1.0);
				} else if(funIndex==2){
					//1 = r + s + t
					return new FConstant(-1.0);
				} else 
					return new FConstant(0.0);
			}
		}
		@Override
		public double value(Variable v) {
			return v.get(sfVarNames.get(funIndex));
		}
		
		public String toString() {
			return sfVarNames.get(funIndex);
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
		if(funID<1 || funID>3) {
			System.out.println("ERROR: funID should be 1,2 or 3.");
			return;
		}
		
		sfVarNames.add("r");
		sfVarNames.add("s");
		sfVarNames.add("t");
		
		//复合函数
		Map<String, FunctionDerivable> fInners = new HashMap<String, FunctionDerivable>();
		List<String> varNamesInner = new LinkedList<String>();
		varNamesInner.add("x");
		varNamesInner.add("y");
		
		final String varName = sfVarNames.get(funIndex);
		fInners.put(varName, new FAbstract(varNamesInner) {	
			public FunctionDerivable derivative(DerivativeIndicator di) {
				if(area < 0.0) {
					FutureEyeException e = new FutureEyeException("SFLinearLocal2D: area < 0.0");
					e.printStackTrace();
					return null;
				}
				if(varName == "r") {
					if(di.get().get("x") != null)
						return new FConstant(b[0] / (2 * area));
					if(di.get().get("y") != null)
						return new FConstant(c[0] / (2 * area));
				} else if(varName == "s") {
					if(di.get().get("x") != null)
						return new FConstant(b[1] / (2 * area));
					if(di.get().get("y") != null)
						return new FConstant(c[1] / (2 * area));
				} else if(varName == "t") {
					if(di.get().get("x") != null)
						return new FConstant(b[2] / (2 * area));
					if(di.get().get("y") != null)
						return new FConstant(c[2] / (2 * area));
				}
				return null;
			}
		});
		
		//使用复合函数构造形函数
		funOuter = new SF123(funIndex);
		this.coef = coef;
		funCompose = FOBasicDerivable.Mult(new FConstant(this.coef), 
				FOBasicDerivable.Compose(funOuter, fInners));
	}
	
	public SFLinearLocal2D(int funID) {
		this.Create(funID, 1.0);
	}
	
	public SFLinearLocal2D(int funID, double coef) {
		this.Create(funID, coef);
	}
	
	@Override
	public FunctionDerivable derivative(DerivativeIndicator di) {
		return funCompose.derivative(di);
	}

	@Override
	public double value(Variable v) {
		return funCompose.value(v);
	}

	@Override
	public void asignElement(Element e) {
		this.e = e;
		//由node改为vertex，因为Element.adjustVerticeToCounterClockwise()结点顺序只调整了vertex
		List<Vertex> vList = e.getVertexList();
		double x1 = vList.get(0).coord(1) , y1 =  vList.get(0).coord(2) ;
		double x2 = vList.get(1).coord(1) , y2 =  vList.get(1).coord(2) ;
		double x3 = vList.get(2).coord(1) , y3 =  vList.get(2).coord(2) ;
		
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

	@Override
	public void setVarNames(List<String> varNames) {
		this.sfVarNames = varNames;
	}

	@Override
	public List<String> varNames() {
		return sfVarNames;
	}
	
	public String toString() {
		String varName = sfVarNames.get(funIndex);
		return "N"+(funIndex+1)+"( "+varName+"(x,y) )="+funOuter.toString();
	}

	
	ShapeFunction sf1d1 = new SFLinearLocal1D(1);
	ShapeFunction sf1d2 = new SFLinearLocal1D(2);
	@Override
	public ShapeFunction restrictTo(int funIndex) {
		if(funIndex == 1) return sf1d1;
		else return sf1d2;
	}
	
}
