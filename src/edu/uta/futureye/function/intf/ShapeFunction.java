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
	 * 当一个形函数与单元关联之后，就可以计算该形函数的函数值和导数值了
	 * JIT编译后的函数：
	 * 选项1:可以有一个成员函数保存该单元
	 * 选项2:删除该函数，从参数直接传入单元坐标
	 * @param e
	 */
	void assignElement(Element e);
	
	/**
	 * 将形函数限制为低一维的形函数
	 * @param funIndex
	 * @return
	 */
	ShapeFunction restrictTo(int funIndex);
	
}
