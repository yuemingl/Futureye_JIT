package edu.uta.futureye.lib.shapefun;

import java.util.HashMap;
import java.util.Map;

import edu.uta.futureye.core.Element;
import edu.uta.futureye.function.AbstractFunction;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.basic.FX;
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
public class SFLinearLocal2DTest  extends AbstractFunction 
							  implements ScalarShapeFunction {
	protected int funIndex;
	private Function funCompose = null;
	private Function funOuter = null;
	protected ObjList<String> innerVarNames = null;
	
	protected Element e = null;
	private double jac = 0.0;
	private double[] x = new double[3];
	private double[] y = new double[3];
	
	private double coef = 1.0;
	
	/**
	 * 构造下列形函数中的一个：
	 * N1 = L1 = r
	 * N2 = L2 = s
	 * N3 = L3 = t = (1 - r - s)
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
		//varNames.add("t"); //DONOT add this name! It's not a free variable
		
		innerVarNames = new ObjList<String>("x","y");
		
		//复合函数
		Map<String, Function> fInners = new HashMap<String, Function>();
		
		for(final String varName : varNames) {
			fInners.put(varName, new AbstractFunction(innerVarNames.toList()) {	
				public Function _d(String var) {
					if(varName.equals("r")) { //r对应三角形高h的负倒数
						if(var.equals("x"))
							return new FC( (y[1]-y[2]) / jac);
						if(var.equals("y"))
							return new FC( (x[2]-x[1]) / jac);
					} else if(varName.equals("s")) { //s对应三角形高h的负倒数
						if(var.equals("x"))
							return new FC( (y[2]-y[0]) / jac);
						if(var.equals("y"))
							return new FC( (x[0]-x[2]) / jac);
					} else
						throw new FutureyeException("\nERROR:\n varName="+varName);
					return null;
				}
				@Override
				public double value(Variable v) {
					throw new FutureyeException("\nERROR:\n Not supported evaluate: "+v);
				}
			});
		}
		
		//Construct shape functions: 
		//r,s are free variables, t=1 - r - s
		if(funIndex == 0)
			funOuter = new FX("r"); //N1=r
		else if(funIndex == 1)
			funOuter = new FX("s"); //N2=s
		else 
			funOuter = FC.c1.S(FX.fr).S(FX.fs); //N3=1-r-s, (N3=t)
		funOuter.setVarNames(varNames);
		
		this.coef = coef;
		funCompose = funOuter.compose(fInners).M(FC.c(this.coef));
		//funCompose.value(new Variable("x",0).set("y",0));
	}
	
	public SFLinearLocal2DTest(int funID) {
		this.Create(funID, 1.0);
	}
	
	public SFLinearLocal2DTest(int funID, double coef) {
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
		VertexList vList = e.vertices();
		x[0] = vList.at(1).coord(1);
		x[1] = vList.at(2).coord(1);
		x[2] = vList.at(3).coord(1);
		y[0] = vList.at(1).coord(2);
		y[1] = vList.at(2).coord(2);
		y[2] = vList.at(3).coord(2);
		jac = (x[0]-x[2])*(y[1]-y[2])-(x[1]-x[2])*(y[0]-y[2]);
	}

	public String toString() {
		return "N"+(funIndex+1)+"(r,s)="+funOuter.toString();
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
