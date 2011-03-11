package edu.uta.futureye.lib.shapefun;

import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

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
import edu.uta.futureye.util.list.VertexList;

public class SFQuadraticLocal1D extends AbstractFunction implements ScalarShapeFunction {
	private int funIndex;
	private Function funCompose = null;
	private Function funOuter = null;
	private List<String> varNames = new LinkedList<String>();
	private ObjList<String> innerVarNames = null;

	private Element e = null;

	/**
	 * 构造下列形函数中的一个：
	 * 
	 *  1--3--2  -->r
	 * -1  0  1
	 * 
	 * N1 = r*(r-1)/2
	 * N2 = (r+1)*r/2
	 * N3 = 1-r*r
	 * @param funID = 1,2,3
	 * 
	 */
	public SFQuadraticLocal1D(int funID) {
		funIndex = funID - 1;
		if(funID<1 || funID>3) {
			System.out.println("ERROR: funID should be 1 ,2 or 3.");
			return;
		}
		
		varNames.add("r");
		innerVarNames = new ObjList<String>("x");
		
		//复合函数
		Map<String, Function> fInners = new HashMap<String, Function>();
		
		/*
		 *  r_x = 2/(x2-x1)  
		 */
		fInners.put("r", new AbstractFunction(innerVarNames.toList()) {	
			public Function d(String var) {
				if(var.contains("x")) {
					VertexList vl = e.vertices();
					if(vl.size() == 2) {
						//TODO ??? 1-0? 0-1?
						double delta = vl.at(2).coord(1)-vl.at(1).coord(1);
						return new FC(2.0/delta);
					} else {
						System.out.println("ERROR: SFLinearLocal1D vl.size()!=2, vl.size()="+vl.size());
					}
				}
				return null;
			}
		});
		
		
		 //* N1 = r*(r-1)/2 = r*(r/2 - 1/2)
		 //* N2 = (r+1)*r/2 = r*(r/2 + 1/2)
		 //* N3 = 1-r*r
		Function fr = new FX("r");
		//使用复合函数构造形函数
		if(funIndex == 0)
			funOuter = FMath.Mult(fr, new FAxpb("r",0.5,-0.5));
		else if(funIndex == 1)
			funOuter = FMath.Mult(fr, new FAxpb("r",0.5,0.5));
		else if(funIndex == 2)
			funOuter = FMath.Minus(new FC(1.0), FMath.Mult(fr, fr));
		
		funCompose = FMath.Compose(funOuter, fInners);
	}
	
	@Override
	public void asignElement(Element e) {
		this.e = e;
	}

	@Override
	public Function d(String varName) {
		return funCompose.d(varName);
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
	public ScalarShapeFunction restrictTo(int funIndex) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public ObjList<String> innerVarNames() {
		return innerVarNames;
	}
}
