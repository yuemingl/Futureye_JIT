package edu.uta.futureye.util.list;

import edu.uta.futureye.core.Node;

/**
 * Node List Class
 * 节点列表类
 * 
 * @author liuyueming
 *
 */
public class NodeList extends ObjList<Node>{
	@Override
	public String toString() {
		return "NodeList"+objs.toString();
	}
}