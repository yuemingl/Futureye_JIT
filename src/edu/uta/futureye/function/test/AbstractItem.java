package edu.uta.futureye.function.test;

import java.util.ArrayList;
import java.util.List;

public abstract class AbstractItem implements Item {
	protected String _name = null;
	protected List<Item> _subItems = new ArrayList<Item>();
	
	@Override
	public String getName() {
		return _name;
	}
	
	@Override
	public void setName(String name) {
		this._name = name;
	}
	
	@Override
	public double getValue(Item ...items) {
		return 0.0;
	}
	
	@Override
	public Item _d(String name) {
		return null;
	}
	
	@Override
	public int symCompairTo(Item item) {
		if(_name == null && item.getName() == null) {
			return 0;
		} else if(_name != null && item.getName() == null) {
			return 1;
		} else if(_name == null && item.getName() != null) {
			return -1;
		} else {
			return item.getName().compareTo(_name);
		}
	}
	
	@Override
	public boolean isDummy() {
		return _name == null;
	}
	
	public String toString() {
		if(_name==null) return "Dummy";
		//if(_subItems.size()>0) return _subItems.toString();
		return _name;
	}
	
	@Override
	public List<Item> getSubItems() {
		return _subItems;
	}
	
	@Override
	public void setSubItems(List<Item> items) {
		_subItems = items;
	}
	
}
