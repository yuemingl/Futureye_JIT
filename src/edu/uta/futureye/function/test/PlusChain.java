package edu.uta.futureye.function.test;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import edu.uta.futureye.util.Constant;

public class PlusChain extends AbstractChain {
	
	@Override
	public double getValue(Item ...items) {
		double rlt = 0.0;
		for(int i=0;i<length();i++) {
			ItemPair pair = getItem(i);
			if(!pair.item.isDummy()) {
				rlt += pair.coef*pair.item.getValue(items);
			} else {
				rlt += pair.coef;
			}
		}
		return rlt;
	}
	
	public void multiCoef(double coef) {
		for(int i=0;i<length();i++) {
			ItemPair pair = getItem(i);
			pair.coef = pair.coef * coef;
		}
	}
	
	@Override
	public void merge(boolean bMergeFunction) {
		List<ItemPair> tmp = new ArrayList<ItemPair>();
		for(int i=0;i<length();i++) {
			ItemPair pair = getItem(i);
			if(bMergeFunction && pair.item instanceof Function) {
				Chain sub = (Chain)((Function)pair.item).getChain().copy();
				sub.merge(bMergeFunction);
				if(sub instanceof PlusChain) {
					((PlusChain) sub).multiCoef(pair.coef);
					tmp.addAll(sub.getAllItem());
				} else {
					//不应该存在
					//Exception e = new Exception("ERROR 1: PlusChain.merge");
					//e.printStackTrace();
					tmp.add(pair);
				}
			} else if(pair.item instanceof PlusChain){
				Chain sub = (Chain)pair.item.copy();
				sub.merge(bMergeFunction);
				((PlusChain) sub).multiCoef(pair.coef);
				tmp.addAll(sub.getAllItem());
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
					if(pair.item.symCompairTo(tmp.get(i).item)==0) {
						pair.coef += tmp.get(i).coef;
						continue;
					}
					tmpMerge.add(pair);
					pair = tmp.get(i);
				} else {
					if(pair.item.getName() == null) {
						pair.coef += tmp.get(i).coef;
					} else if(pair.item.getName().equals(tmp.get(i).item.getName())) {
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
		List<ItemPair> tmp = new ArrayList<ItemPair>();
		for(int i=0;i<length();i++) {
			ItemPair pair = getItem(i);
			Item dItem = pair.item._d(name);
			if(dItem != null) {
				tmp.add(new ItemPair(pair.coef,dItem));
			}
		}
		PlusChain pc = new PlusChain();
		pc.addAllItem(tmp);
		pc.merge(false);
		return pc;
	}

	public Item copy() {
		Chain c = new PlusChain();
		for(ItemPair pair : this._items) {
			ItemPair newPair = new ItemPair(pair.coef,pair.item.copy());
			c.addItem(newPair);
		}
		return c;
	}	
	
	public String toString() {
		String rlt = "  ";
		String operator = "+";
		
		for(int i=0;i<length();i++) {
			ItemPair pair = getItem(i);
			if(pair.item instanceof Chain) {
				rlt += pair.coef+"*"+pair.item.toString()+" "+operator+" ";
			} else {
				if(!pair.item.isDummy()) {
					if(Math.abs(pair.coef-1.0) < Constant.eps) {
						rlt += pair.item.toString()+" "+operator+" ";
					} else {
						if(pair.item instanceof Function)
							rlt += pair.coef+"*("+((Function)pair.item).toString()+") "+operator+" ";
						else
							rlt += pair.coef+"*"+pair.item.toString()+" "+operator+" ";
					}
				} else {
					rlt += pair.coef+" "+operator+" ";;
				}
			}
		}
		return rlt.substring(0,rlt.length()-2);	
	}
	

	
}
