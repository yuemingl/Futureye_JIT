package edu.uta.futureye.lib.shapefun;

import edu.uta.futureye.core.Element;
import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.intf.ScalarShapeFunction;
import edu.uta.futureye.function.intf.ShapeFunction;
import edu.uta.futureye.util.container.ObjList;

public class SF0 extends FC implements ScalarShapeFunction {
	protected ObjList<String> innerVarNames = new ObjList<String>();

	public SF0(String ...varNames) {
		for(int i=1;i<varNames.length;i++) {
			this.innerVarNames.add(varNames[i]);
		}
	}
	public SF0(ObjList<String> varNames) {
		for(int i=1;i<=varNames.size();i++) {
			this.innerVarNames.add(varNames.at(i));
		}
	}
	
	@Override
	public void asignElement(Element e) {
	}

	@Override
	public ObjList<String> innerVarNames() {
		return innerVarNames;
	}

	@Override
	public ShapeFunction restrictTo(int funIndex) {
		return null;
	}
}
