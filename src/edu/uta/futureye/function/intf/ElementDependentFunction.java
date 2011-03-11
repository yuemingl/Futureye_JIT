package edu.uta.futureye.function.intf;

import edu.uta.futureye.core.Element;

/**
 * 与单元相关的函数，例如：一个函数的方向导数与单元的边界法方向有关，
 * 因此，需要单元信息
 * 
 * @author liuyueming
 *
 */
public interface ElementDependentFunction extends Function {
	void setElement(Element e);
}
