package edu.uta.futureye.function.test;

public class ISin extends AbstractItem {
	public ISin() {
		_subItems.add(new Variable("x"));
		setName("sin(x)");
	}

	public ISin(String varName) {
		_subItems.add(new Variable("x"));
		setName("sin("+varName+")");
	}
	
	public ISin(Variable var) {
		_subItems.add(new Variable("x"));
		setName("sin("+var.getName()+")");
	}
	
	@Override
	public double getValue(Item ...items) {
		if(items == null || items.length==0) {
			Exception e = new Exception("ISin getValue");
			e.printStackTrace();
		} else {
			for(int i=0;i<items.length;i++) {
				for(int j=0;j<_subItems.size();j++) {
					if(_subItems.get(j).getName().equals(items[i].getName())) {
						return Math.sin(
							_subItems.get(j).getValue(items)
						);
					}
				}
			}
 		}
		return 0.0;
	}
	
	@Override
	public Item _d(String name) {
		for(int j=0;j<_subItems.size();j++) {
			if(_subItems.get(j).getName().equals(name)) {
				Function fc = new FCos(name);
				fc.setSubItems(this._subItems);
				return fc;
			} else if(_subItems.get(j) instanceof Function) {
				Function fc = new FCos(name);
				fc.getChain().getItem(0).item.setSubItems(this._subItems);
				//fc.setSubItems(this._subItems);
				Function fd = (Function)_subItems.get(j)._d(name);
				if(fd != null) {
					return FOperator.Multi(fc,fd);
				}
			}
		}
		return null;
	}
	
	@Override
	public Item copy() {
		Item i = new ISin();
		i.setName(this.getName());
		i.setSubItems(this.getSubItems());
		return i;
	}
	
	public String toString() {
		if(_subItems.size()>0) {
			String s = _subItems.toString();
			return "sin( "+s.substring(1,s.length()-1)+" )";
		}
		return _name;
	}
}
