package edu.uta.futureye.algebra.intf;

import java.util.Map;

public interface SparseMatrix extends Matrix,Iterable<MatrixEntry> {
	/**
	 * 
	 * @return
	 */
	int getNonZeroNumber();
	
	
	/**
	 * get all non-zero element, instead of iterator
	 * <p>
	 * 获取所有非零元素，不使用迭代子
	 * @return
	 */
	Map<Integer,Map<Integer,Double>> getAll();
	
	
	/**
	 *TODO ??? setPart()
	 * 
	 * M(nRowBase + row, nColBase + col) = values in map
	 * <p>
	 * 将参数map中的所有元素值赋值到本矩阵，其中map中所有的行号都加上nRowBase，
	 * 所有的列号都加上nColBase
	 * @param nRowBase
	 * @param nColBase
	 * @param dataMap
	 */
	void setAll(int nRowBase, int nColBase, Map<Integer,Map<Integer,Double>> dataMap);
	
	
	/**
	 * Clear data and dimension information
	 * <p>
	 * 清空矩阵中的所有数据和维度信息
	 */
	void clearAll();
	
	/**
	 * Clear data only, keep matrix dimension information
	 * <p>
	 * 只清空矩阵中的数据，保留维度信息
	 */
	void clearData();
	
	/**
	 * An overriding method can also return a subtype of the type returned by the overridden method. 
	 * This is called a covariant return type.
	 */
	SparseMatrix trans();
	SparseMatrix copy();
	SparseMatrix setName(String name);
}
