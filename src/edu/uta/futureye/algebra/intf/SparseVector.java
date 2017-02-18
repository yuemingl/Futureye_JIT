/**
 * Copyright (c) 2010, nkliuyueming@gmail.com. All rights reserved.
 * 
 * 
 */
package edu.uta.futureye.algebra.intf;

import java.util.Map;

public interface SparseVector extends Vector,Iterable<VectorEntry> {
	/**
	 * 
	 * @return
	 */
	int getNonZeroNumber();
	
	/**
	 * 
	 * @return
	 */
	Map<Integer,Double> getAll();
	
	/**
	 * 
	 * @param dataMap
	 * @return
	 */
	SparseVector setAll(int nBase, Map<Integer,Double> dataMap);
	
	/**
	 * Clear all values but keep dimension of vector
	 * 清空向量的所有元素，但保留维度信息
	 */
	void clearData();
	
	/**
	 * Clear all values and set dimension of vector to zero.
	 * <p>
	 * 清空向量的所有元素并且维度置零
	 */
	void clearAll();
	
	/**
	 * An overriding method can also return a subtype of the type returned by the overridden method. 
	 * This is called a covariant return type.
	 */
	SparseVector setName(String name);
	SparseVector copy();
	SparseVector set(Vector y);
	SparseVector set(double a, Vector y);
	SparseVector add(Vector y);
	SparseVector add(double a, Vector y);
	SparseVector scale(double a);
	SparseVector ax(double a);
	SparseVector axpy(double a, Vector y);
	SparseVector axMuly(double a, Vector y);
	SparseVector axDivy(double a, Vector y);
	SparseVector shift(double dv);
}
