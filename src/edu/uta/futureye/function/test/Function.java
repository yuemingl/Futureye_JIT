package edu.uta.futureye.function.test;

public interface Function extends Item {

	/**
	 * 创建函数表达式
	 */
	public void createChain();
	
	/**
	 *  获取表达式
	 * @return
	 */
	public Chain getChain();
	
	/**
	 * 设置表达式
	 */
	public void setChain(Chain chain);
	
	/**
	 * 复合函数
	 * @param f
	 * @return
	 */
	public Function compose(ComposePair ...pairs);
	
	/**
	 * 展开多项式
	 * @return
	 */
	public Function expand();	
	
	/**
	 * 函数求值（多自变量）
	 * @param v
	 * @return
	 */
	public double getValue(Variable ...v);
}
