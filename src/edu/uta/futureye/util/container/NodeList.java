/**
 * Copyright (c) 2010, nkliuyueming@gmail.com. All rights reserved.
 * 
 * 
 */
package edu.uta.futureye.util.container;

import edu.uta.futureye.core.Node;

/**
 * <P>Node List Container</P>
 * <P>节点列表容器</P>
 * 
 * <B>Notes:</B>
 * <P>Node index starts from 1.</P>
 * <P><tt>null</tt> element is not allowed.</P> 
 * <P>Auto size increment is supported.</P>

 * @author liuyueming
 */
public class NodeList extends ObjList<Node>{
	@Override
	public String toString() {
		return "NodeList"+objs.toString();
	}
}