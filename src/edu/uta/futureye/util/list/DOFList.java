package edu.uta.futureye.util.list;

import edu.uta.futureye.core.DOF;

/**
 * Degree Of Freedom Class
 * 自由度列表类
 * 
 * @author liuyueming
 *
 */
public class DOFList extends ObjList<DOF> {
	@Override
	public String toString() {
		return "DOFList"+objs.toString();
	}
}
