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
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.function.intf.ScalarShapeFunction;
import edu.uta.futureye.function.operator.FMath;
import edu.uta.futureye.util.list.ObjList;
import edu.uta.futureye.util.list.VertexList;

public class SFLinearLocal1D extends AbstractFunction  implements ScalarShapeFunction {
	private int funIndex;
	private Function funCompose = null;
	private Function funOuter = null;
	private List<String> varNames = new LinkedList<String>();
	private ObjList<String> innerVarNames = null;

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
		innerVarNames = new ObjList<String>("x");
		
		//复合函数
		Map<String, Function> fInners = new HashMap<String, Function>();
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
		fInners.put("r", new AbstractFunction(varNamesInner) {	
			public Function d(String var) {
				if(var.equals("x")) {
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
		
		//使用复合函数构造形函数
		if(funIndex == 0)
			funOuter = new FAxpb("r",-0.5,0.5);
		else
			funOuter = new FAxpb("r",0.5,0.5);
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
