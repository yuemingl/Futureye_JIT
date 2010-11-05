package edu.uta.futureye.util;

import java.util.ArrayList;
import java.util.List;

import edu.uta.futureye.core.Element;

/**
 * Element List Class
 * 单元列表类
 * @author liuyueming
 * 
 */
public class ElementList {
	protected List<Element> elements = new ArrayList<Element>();
	
	/**
	 * @param index start from 1,2,3...
	 * @return
	 */	
	public Element at(int index) {
		if(index < 1)
			System.out.println("ERROR: ElementList index="+index);
		return elements.get(index-1);
	}

	public void add(Element e) {
		this.elements.add(e);
	}
	
	public int size() {
		return elements.size();
	}

	public void clear() {
		elements.clear();
	}
	
	public String toString() {
		return elements.toString();
	}
	
	public Element remove(int index) {
		return elements.remove(index-1);
	}
	
	public boolean remove(Element e) {
		return elements.remove(e);
	}
	
	public void addAll(ElementList list) {
		for(int i=1;i<=list.size();i++)
			elements.add(list.at(i));
	}
	
}
