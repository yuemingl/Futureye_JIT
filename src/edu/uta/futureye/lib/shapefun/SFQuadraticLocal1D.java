package edu.uta.futureye.lib.shapefun;

import java.util.HashMap;
import java.util.Map;

import edu.uta.futureye.core.Element;
import edu.uta.futureye.function.AbstractMathFunc;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.FAxpb;
import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.basic.FX;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.function.intf.ScalarShapeFunction;
import edu.uta.futureye.util.FutureyeException;
import edu.uta.futureye.util.container.ObjList;
import edu.uta.futureye.util.container.VertexList;

public class SFQuadraticLocal1D extends AbstractMathFunc implements ScalarShapeFunction {
	private int funIndex;
	private MathFunc funCompose = null;
	private MathFunc funOuter = null;
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
			FutureyeException ex = new FutureyeException("ERROR: funID should be 1 ,2 or 3.");
			ex.printStackTrace();
			System.exit(-1);			
		}
		
		varNames.add("r");
		innerVarNames = new ObjList<String>("x");
		
		//复合函数
		Map<String, MathFunc> fInners = new HashMap<String, MathFunc>();
		
		/*
		 *  r_x = 2/(x2-x1)  
		 */
		fInners.put("r", new AbstractMathFunc(innerVarNames.toList()) {	
			public MathFunc diff(String var) {
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

			@Override
			public double apply(Variable v) {
				// TODO Auto-generated method stub
				return 0;
			}
		});
		
		
		 //* N1 = r*(r-1)/2 = r*(r/2 - 1/2)
		 //* N2 = (r+1)*r/2 = r*(r/2 + 1/2)
		 //* N3 = 1-r*r
		MathFunc fr = new FX("r");
		//使用复合函数构造形函数
		if(funIndex == 0)
			funOuter = fr.M(new FAxpb("r",0.5,-0.5));
		else if(funIndex == 1)
			funOuter = fr.M(new FAxpb("r",0.5,0.5));
		else if(funIndex == 2)
			funOuter = FC.C1.S(fr.M(fr));
		
		funCompose = funOuter.compose(fInners);
	}
	
	@Override
	public void assignElement(Element e) {
		this.e = e;
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
