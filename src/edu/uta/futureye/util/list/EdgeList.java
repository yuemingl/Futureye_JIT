package edu.uta.futureye.util.list;

import edu.uta.futureye.core.Edge;

/**
 * Edge List Class
 * （全局）边类列表
 * 
 * @author liuyueming
 *
 */
public class EdgeList extends ObjList<Edge> {
	@Override
	public String toString() {
		return "EdgeList"+objs.toString();
	}
}

