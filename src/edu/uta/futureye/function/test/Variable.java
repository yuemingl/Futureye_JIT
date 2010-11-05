package edu.uta.futureye.function.test;

public class Variable extends AbstractItem {
	double _value = 0.0;
	
	public Variable(String name) {
		this._name = name;
	}
	
	public Variable(String name, double value) {
		this._name = name;
		this._value = value;
	}
	
	@Override
	public double getValue(Item ...items) {
		if(items == null || items.length==0)
			return _value;
		else {
			for(int i=0;i<items.length;i++) {
				if(_name.equals(items[i].getName())) {
					return ((Variable)items[i]).getValue();
				}
			}
 		}
		return _value;
	}
	
	@Override
	public Item _d(String name) {
		if(name.equals(this._name)) {
			return new GhostItem();
		}
		return null;
	}

	@Override
	public Item copy() {
		Variable v = new Variable(this._name,this._value);
		return v;
	}
}
