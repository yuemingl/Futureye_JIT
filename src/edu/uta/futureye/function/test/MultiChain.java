package edu.uta.futureye.function.test;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import edu.uta.futureye.util.Constant;

public class MultiChain extends AbstractChain {
	
	@Override
	public double getValue(Item ...items) {
		double rlt = 1.0;
		for(int i=0;i<length();i++) {
			ItemPair pair = getItem(i);
			//pair.coef 指数阶
			if(!pair.item.isDummy()) {
				rlt *= Math.pow(pair.item.getValue(items), pair.coef);
			} else {
				rlt *= pair.coef;
			}
		}
		return rlt;
	}
	
	public void multiPower(double power) {
		for(int i=0;i<length();i++) {
			ItemPair pair = getItem(i);
			pair.coef = pair.coef * power;
		}
	}
	
	@Override
	public void merge(boolean bMergeFunction) {
		if(length()<=1) return;
		List<ItemPair> tmp = new ArrayList<ItemPair>();
		for(int i=0;i<length();i++) {
			ItemPair pair = getItem(i);
			if(bMergeFunction && pair.item instanceof Function) {
				Chain sub = ((Function)pair.item).getChain();
				sub.merge(bMergeFunction);
				if(sub instanceof MultiChain) {
					tmp.addAll(sub.getAllItem());
				} else {
					tmp.add(pair);
				}
			} else if(pair.item instanceof MultiChain) {
				tmp.addAll(((MultiChain)pair.item).getAllItem());
			} else {
				tmp.add(pair);
			}
		}
		Collections.sort(tmp,new Comparator<ItemPair>(){
			@Override
			public int compare(ItemPair o1, ItemPair o2) {
				if(o1.item instanceof Chain || o1.item instanceof Function)
					return -o1.item.symCompairTo(o2.item);
				else
					return o2.item.symCompairTo(o1.item);
			}
		});
		this._items.clear();
		List<ItemPair> tmpMerge = new ArrayList<ItemPair>();
		if(tmp.size() == 1) {
			this._items = tmp;
			return;
		} else {
			ItemPair pair = tmp.get(0);
			for(int i=1;i<tmp.size();i++) {
				if(pair.item instanceof Chain) {
					tmpMerge.add(pair);
					pair = tmp.get(i);
				} else {
					if(pair.item.getName() == null) {
						//系数部分相乘
						pair.coef *= tmp.get(i).coef;
					} else if(pair.item.getName().equals(tmp.get(i).item.getName())) {
						//指数部分相加
						pair.coef += tmp.get(i).coef;
					} else {
						tmpMerge.add(pair);
						pair = tmp.get(i);
					}
				}
			}
			tmpMerge.add(pair);
		}
		this._items = tmpMerge;
	}
	
	@Override
	public Item _d(String name) {
		PlusChain pc = new PlusChain();
		for(int i=0;i<length();i++) {
			MultiChain mc = new MultiChain();
			ItemPair pair = getItem(i);
			double coef = pair.coef;
			Item dItem = null;
			if(Math.abs(pair.coef-1.0)<Constant.eps) {
				dItem = pair.item._d(name);
				if(dItem != null)
					mc.addItem(new ItemPair(1.0,dItem));
				else
					mc.addItem(new ItemPair(0.0));
			} else {
				dItem = pair.item._d(name);
				if(dItem != null) {
					mc.addItem(new ItemPair(pair.coef-1,pair.item));
					mc.addItem(new ItemPair(1.0,dItem));
				} else {
					mc.addItem(new ItemPair(0.0));
				}
			}
			for(int j=0;j<length();j++) {
				ItemPair pair2 = getItem(j);
				if(i != j) {
					mc.addItem(pair2);
				}
			}
			//系数
			mc.addItem(new ItemPair(coef));
			mc.merge(false);
			pc.addItem(new ItemPair(1.0,mc));
		}
		pc.merge(false);
		return pc;
	}
	
	public Item copy() {
		Chain c = new MultiChain();
		for(ItemPair pair : this._items) {
			ItemPair newPair = new ItemPair(pair.coef,pair.item.copy());
			c.addItem(newPair);
		}
		return c;
	}	
	
	public String toString() {
		String rlt = "  ";
		for(int i=0;i<length();i++) {
			ItemPair pair = getItem(i);
			if(pair.item instanceof Chain) {
				if(Math.abs(pair.coef-1.0) < Constant.eps)
					rlt += " ("+pair.item.toString()+") * ";
				else
					rlt += " ("+pair.item.toString()+")^"+pair.coef+" * ";
			} else {
				if(!pair.item.isDummy()) {
					if(Math.abs(pair.coef-1.0) < Constant.eps)
						rlt += " ("+pair.item.toString()+") * ";
					else
						rlt += " ("+pair.item.toString()+")^"+pair.coef+" * ";
				} else {
					rlt += pair.coef + " * ";
				}
			}
		}
		return rlt.substring(0,rlt.length()-2);	
	}
	
	public ItemPair reduceToItemPair() {
		List<ItemPair> items = new ArrayList<ItemPair>();
		double coef = 1.0;
		for(int i=0;i<length();i++) {
			ItemPair pair = getItem(i);
			if(pair.item.isDummy()) {
				coef *= pair.coef;
			} else {
				items.add(pair);
			}
		}
		if(items.size() == 0) {
			return new ItemPair(coef);
		} else if(items.size() == 1) {
			ItemPair pair = items.get(0);
			if(Math.abs(pair.coef - 1.0) < Constant.eps)
				return new ItemPair(coef,pair.item);
			else {
				this._items = items;
				return new ItemPair(coef,this);
			}
		} else {
			this._items = items;
			return new ItemPair(coef,this);
		}
	}
	
}
