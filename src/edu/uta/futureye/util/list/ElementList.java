package edu.uta.futureye.util.list;

import edu.uta.futureye.core.Element;

/**
 * Element List Class
 * 单元列表类
 * 
 * @author liuyueming
 * 
 */
public class ElementList extends ObjList<Element>{
	@Override
	public String toString() {
		return "ElementList"+objs.toString();
	}
}
