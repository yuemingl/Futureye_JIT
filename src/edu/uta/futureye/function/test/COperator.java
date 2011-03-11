package edu.uta.futureye.function.test;

import java.util.List;

import edu.uta.futureye.util.Constant;

public class COperator {
	
	/**
	 * item * pChain
	 * @param pChain
	 * @param item
	 * @return
	 */
	public static PlusChain PlusChainMultiItem(PlusChain pChain, Item item) {
		if(item instanceof MultiChain)
			return PlusChainMultiMultiChain(pChain,(MultiChain)item);
		else if(item instanceof PlusChain)
			return PlusChainMultiPlusChain(pChain,(PlusChain)item);
		else {
			PlusChain pc = new PlusChain();
			for(ItemPair pair : pChain.getAllItem()) {
				MultiChain mc = new MultiChain();
				if(Math.abs(pair.coef-1.0) > Constant.eps)
					mc.addItem(new ItemPair(pair.coef));
				mc.addItem(new ItemPair(1.0,pair.item));
				mc.addItem(new ItemPair(1.0,item));
				pc.addItem(mc.reduceToItemPair());
			}
			return pc;
		}
	}
	
	/**
	 * pChain * mChain
	 * @param pChain
	 * @param mChain
	 * @return
	 */
	public static PlusChain PlusChainMultiMultiChain(PlusChain pChain, MultiChain mChain) {
		PlusChain pc = new PlusChain();
		for(ItemPair pair : pChain.getAllItem()) {
			MultiChain mc = new MultiChain();
			if(Math.abs(pair.coef-1.0) > Constant.eps)
				mc.addItem(new ItemPair(pair.coef));
			mc.addItem(new ItemPair(1.0,pair.item));
			mc.addAllItem(mChain.getAllItem());
			pc.addItem(mc.reduceToItemPair());
		}
		return pc;
	}
	
	/**
	 * pChain1 * pChain2
	 * @param pChain1
	 * @param pChain2
	 * @return
	 */
	public static PlusChain PlusChainMultiPlusChain(PlusChain pChain1, PlusChain pChain2) {
		PlusChain pc = new PlusChain();
		for(ItemPair pair1 : pChain1.getAllItem()) {
			for(ItemPair pair2 : pChain2.getAllItem()) {
				MultiChain mc = new MultiChain();
				double coef = pair1.coef * pair2.coef;
				if(Math.abs(coef-1.0) > Constant.eps)
					mc.addItem(new ItemPair(coef));
				mc.addItem(new ItemPair(1.0,pair1.item));
				mc.addItem(new ItemPair(1.0,pair2.item));
				pc.addItem(mc.reduceToItemPair());
			}
		}
		return pc;
	}
	
	
	public static Chain ReduceChain(Chain chain) {
		chain.merge(true);
		if(chain instanceof PlusChain) {
			for(int i=0;i<chain.length();i++) {
				ItemPair pair = chain.getItem(i);
				if(Math.abs(pair.coef)<Constant.eps)
					pair.item = new GhostItem();
			}
		} else if(chain instanceof MultiChain) {
			for(int i=0;i<chain.length();i++) {
				ItemPair pair = chain.getItem(i);
				if(Math.abs(pair.coef)<Constant.eps && !(pair.item instanceof GhostItem)) {
					pair.coef = 1.0;
					pair.item = new GhostItem();
				}
			}
			chain.merge(false);
		}		
		if(chain.length() == 1) {
			if(chain instanceof PlusChain) {
				ItemPair pair = chain.getItem(0);
				Item item = pair.item;
				if(pair.item instanceof Function) {
					Chain sub = ((Function)pair.item).getChain();
					item = sub;
				}
				if(item instanceof PlusChain) {
					PlusChain pc = (PlusChain)((PlusChain)item).copy();
					pc.multiCoef(pair.coef);
					return ReduceChain(pc);
				} else if(item instanceof MultiChain) {
					MultiChain mc = (MultiChain)((MultiChain)item).copy();
					mc.addItem(new ItemPair(pair.coef));
					return ReduceChain(mc);
				}
			} else if(chain instanceof MultiChain) {
				ItemPair pair = chain.getItem(0);
				Item item = pair.item;
				if(pair.item instanceof Function) {
					Chain sub = ((Function)pair.item).getChain();
					item = sub;
				}
				if(item instanceof PlusChain) {
					if(Math.abs(pair.coef-1.0) < Constant.eps) {
						return ReduceChain((Chain)item);
					}
				} else if(item instanceof MultiChain) {
					if(Math.abs(pair.coef-1.0) < Constant.eps)
						return ReduceChain((Chain)item);
				}
			}
		}

		return chain;
	}
	
	
	/**
	 * 
	 * @param chain
	 * @return
	 */
	public static Chain ExpandChain(Chain chain) {
		chain.merge(true);
		chain = ReduceChain(chain);
		if(chain instanceof PlusChain) {
			PlusChain tmp = (PlusChain)chain.copy();
			List<ItemPair> list = tmp.getAllItem();
			for(int i=0;i<list.size();i++) {
				ItemPair pair = list.get(i);
				if(pair.item instanceof Chain)
					pair.item = ReduceChain((Chain)pair.item);	
				if(pair.item instanceof MultiChain) {
					pair.item = ExpandChain((Chain)pair.item);//Update Item
					if(pair.item instanceof MultiChain) {
						ItemPair reduced = ((MultiChain)pair.item).reduceToItemPair();
						pair.coef *= reduced.coef;
						pair.item = reduced.item;
					} else if(pair.item instanceof PlusChain) {
						pair.item = ExpandChain((Chain)pair.item);
					} else {
						//do nothing
					}
				} else if(pair.item instanceof PlusChain) {
					Exception e = new Exception("ERROR 1: CheckChain");
					e.printStackTrace();
				} else { //item
					//do nothing
				}
			}
			chain = tmp;
		} else if(chain instanceof MultiChain) {
			MultiChain tmp = (MultiChain)chain.copy();
			List<ItemPair> list = tmp.getAllItem();
			for(int i=0;i<list.size();i++) {
				ItemPair pair = list.get(i);
				if(pair.item instanceof Chain)
					pair.item = ReduceChain((Chain)pair.item);
				if(pair.item instanceof PlusChain) {
					if(Math.abs(pair.coef-1.0) < Constant.eps) {
						Item tmp2 = tmp;
						list.remove(i);
						if(list.size() == 1) {
							if(Math.abs(list.get(0).coef-1.0) < Constant.eps) {
								tmp2 = list.get(0).item;
							} else {
								//do nothing
							}
						}
						PlusChain pc = null;
						if(tmp.length()>0) {
							pc = PlusChainMultiItem((PlusChain) pair.item,tmp2);
						} else {
							pc = (PlusChain)pair.item;
							pc.multiCoef(pair.coef);
						}
						return ExpandChain(pc);
					} else if(pair.coef > 0 && Math.abs(Math.ceil(pair.coef)-pair.coef) < Constant.eps) {
						//整数次幂转换成连乘形式
						pair.coef = 1.0;
						for(int j=1;j<=(int)pair.coef;j++) {
							tmp.addItem(new ItemPair(1.0,pair.item));
						}
						return ExpandChain(tmp);
					} else {
						//do nothing
					}
				} else{
					//do nothing
				}
			}
		} else {
			Exception e = new Exception("ERROR: Unsupported Chain Type");
			e.printStackTrace();
		}
		chain.merge(true);
		chain = ReduceChain(chain);
		return chain;
	}
}
