package edu.uta.futureye.function.intf;

import edu.uta.futureye.core.Element;
import edu.uta.futureye.util.container.ObjList;

public interface ShapeFunction {
	/**
	 * 当形函数是局部坐标的复合函数时，返回物理坐标变量名称
	 * @return
	 */
	ObjList<String> innerVarNames();
	
	/**
	 * 关联形函数和单元
	 * @param e
	 */
	void asignElement(Element e);
	
	/**
	 * 将形函数限制为低一维的形函数
	 * @param funIndex
	 * @return
	 */
	ShapeFunction restrictTo(int funIndex);
	
}
