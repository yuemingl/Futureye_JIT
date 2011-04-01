package edu.uta.futureye.function.test;

import java.util.List;

public interface Item {
	public void setName(String name);
	public String getName();
	
	/**
	 * 从符号上判断两个Item对象代表的“函数项”的关系
	 * @param item
	 * @return 0相等 ，-1小于， 1大于
	 */
	public int symCompairTo(Item item);
	
	/**
	 * 是否亚元，一般情况下，每个item都有系数，
	 * 如果仅表示系数时，可以使用一个Dummy Item
	 * @return
	 */
	public boolean isDummy();
	
	/**
	 * 复制
	 * @return
	 */
	public Item copy();
	
	/**
	 * 求值
	 * @param items
	 * @return
	 */
	public double getValue(Item ...items);
	
	/**
	 * 关于自变量name求导一次
	 * @param name
	 * @return
	 */
	public Item _d(String name);
	
	/**
	 * 
	 * @return
	 */
	public List<Item> getSubItems();
	public void setSubItems(List<Item> items);
}
