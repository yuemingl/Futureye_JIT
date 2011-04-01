package edu.uta.futureye.function.test;

import java.util.ArrayList;
import java.util.List;

import edu.uta.futureye.util.Constant;

public abstract class AbstractChain extends AbstractItem implements Chain {
	protected List<ItemPair> _items = new ArrayList<ItemPair>();
	
	@Override
	public int length() {
		return _items.size();
	}
	
	@Override
	public void addItem(ItemPair pair) {
		_items.add(pair);
	}

	@Override
	public void addAllItem(List<ItemPair> list) {
		_items.addAll(list);
	}

	@Override
	public void setItem(int index,ItemPair pair) {
		_items.set(index, pair);
	}

	@Override
	public List<ItemPair> getAllItem() {
		return _items;
	}

	@Override
	public ItemPair getItem(int index) {
		return _items.get(index);
	}
	
	public String toString() {
		return _items.toString();
	}

	@Override
	public boolean isDummy() {
		return false;
	}
	
	@Override
	public void clear() {
		_items.clear();
	}

	@Override
	public int symCompairTo(Item item) {
		if(item instanceof Chain) {
			Chain c = (Chain)item;
			this.merge(false);
			c.merge(false);
			int minLength = Math.min(this.length(), c.length());
			//默认假设每个Chain中的Item都是从大到小排列的
			for(int i=0;i<minLength;i++) {
				int cmpRlt = this.getItem(i).item.symCompairTo(c.getItem(i).item);
				if(cmpRlt == 0) {
					if(Math.abs(this.getItem(i).coef-c.getItem(i).coef) < Constant.eps) {
						continue;
					} else {
						return this.getItem(i).coef-c.getItem(i).coef>0?1:-1;
					}
				} else {
					return cmpRlt;
				}
			}
			if(this.length() > c.length())
				return 1;
			else if(this.length() < c.length())
				return -1;
			else {
				return 0;
			}
		} else {
			return 1;
		}
	}
}
