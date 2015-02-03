package edu.uta.futureye.lib.shapefun;

import edu.uta.futureye.core.Element;
import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.intf.ScalarShapeFunction;
import edu.uta.futureye.function.intf.ShapeFunction;
import edu.uta.futureye.util.container.ObjList;

/**
 * Shape function: constant zero, as zero component of vector valued shape function
 * <p>
 * 实现形函数接口，返回常值0，作为向量值形函数的0分量
 * 
 * @author liuyueming
 *
 */
public class SFConstant0 extends FC implements ScalarShapeFunction {
	protected ObjList<String> innerVarNames = new ObjList<String>();

	public SFConstant0(String ...varNames) {
		super(0.0);
		for(int i=1;i<varNames.length;i++) {
			this.innerVarNames.add(varNames[i]);
		}
	}
	public SFConstant0(ObjList<String> varNames) {
		super(0.0);		
		for(int i=1;i<=varNames.size();i++) {
			this.innerVarNames.add(varNames.at(i));
		}
	}
	
	@Override
	public void assignElement(Element e) {
	}

	@Override
	public ObjList<String> innerVarNames() {
		return innerVarNames;
	}

	@Override
	public ShapeFunction restrictTo(int funIndex) {
		return null;
	}
	
	@Override
	public boolean isConstant() {
		return true;
	}
}
