package edu.uta.futureye.function.test;

import java.util.List;

public interface Chain extends Item{
	public int length();
	
	public void addItem(ItemPair pair);
	public void addAllItem(List<ItemPair> list);
	public void setItem(int index,ItemPair pair);
	public ItemPair getItem(int index);
	public List<ItemPair> getAllItem();
	public void clear();
	
	/**
	 * 合并Chain中的同类项
	 * @param bMergeFunction 控制是否合并Chain中的Function对象
	 */
	public void merge(boolean bMergeFunction);
	
}
