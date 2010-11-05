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
import edu.uta.futureye.function.basic.FAxpb;
import edu.uta.futureye.function.basic.FConstant;
import edu.uta.futureye.function.intf.FunctionDerivable;
import edu.uta.futureye.function.intf.ShapeFunction;
import edu.uta.futureye.function.operator.FOBasicDerivable;

public class SFLinearLocal1D implements ShapeFunction {
	private int funIndex;
	private FunctionDerivable funCompose = null;
	private FunctionDerivable funOuter = null;
	private List<String> varNames = new LinkedList<String>();

	private Element e = null;

	/**
	 * 构造下列形函数中的一个：
	 * 
	 *  1-----2  -->r
	 * -1  0  1
	 * 
	 * N1 = (1-r)/2
	 * N2 = (1+r)/2
	 * @param funID = 1,2
	 * 
	 */
	public SFLinearLocal1D(int funID) {
		funIndex = funID - 1;
		if(funID<1 || funID>2) {
			System.out.println("ERROR: funID should be 1 or 2.");
			return;
		}
		
		varNames.add("r");
		
		//复合函数
		Map<String, FunctionDerivable> fInners = new HashMap<String, FunctionDerivable>();
		List<String> varNamesInner = new LinkedList<String>();
		varNamesInner.add("x");
		
		/*
		 *  x = x1*N1 + x2*N2
		 *    = x1*(1-r)/2 + x2*(1+r)/2
		 *    = [ x1+x2 + (x2-x1)*r ]/2
		 *  =>
		 *  r = [2*x - (x1+x2)]/(x2-x1) 
		 *  r_x = 2/(x2-x1)  
		 */
		fInners.put("r", new FAbstract(varNamesInner) {	
			public FunctionDerivable derivative(DerivativeIndicator di) {
				
				if(di.get().get("x") != null) {
					List<Vertex> vl = e.getVertexList();
					if(vl.size() == 2) {
						//TODO ??? 1-0? 0-1?
						double delta = vl.get(1).coord(1)-vl.get(0).coord(1);
						return new FConstant(2.0/delta);
					} else {
						System.out.println("ERROR: SFLinearLocal1D vl.size()!=2, vl.size()="+vl.size());
					}
				}
				return null;
			}
		});
		
		//使用复合函数构造形函数
		if(funIndex == 0)
			funOuter = new FAxpb("r",-0.5,0.5);
		else
			funOuter = new FAxpb("r",0.5,0.5);
		funCompose = FOBasicDerivable.Compose(funOuter, fInners);
	}
	
	@Override
	public void asignElement(Element e) {
		this.e = e;
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

	@Override
	public ShapeFunction restrictTo(int funIndex) {
		// TODO Auto-generated method stub
		return null;
	}

}
