package edu.uta.futureye.algebra.intf;

import java.util.Map;

/**
 * General matrix interface
 * 一般矩阵接口，不适用于求解代数方程组（求解代数方程组的矩阵接口参见AlgebraMatrix）
 * 用途：
 * 1. 用来保存刚度矩阵合成步骤中的局部和全局矩阵
 * 2. 打印输出矩阵
 * 
 * @author liuyueming
 *
 */
public interface Matrix {
	/**
	 * 零元素阈值
	 */
	public static double zeroEps = 1e-10;
	
	/**
	 * 设置矩阵行数
	 * @param nRowDim TODO
	 */
	void setRowDim(int nRowDim);
	
	/**
	 * 获取矩阵行数
	 * @return TODO
	 */
	int getRowDim();
	
	/**
	 * 设置矩阵列数
	 * @param nColDim TODO
	 */
	void setColDim(int nColDim);
	
	/**
	 * 获取矩阵列数
	 * @return
	 */
	int getColDim();
	
	/**
	 * 设置row行col列元素的值
	 * @param row
	 * @param col
	 * @param value
	 */
	void set(int row, int col,double value);
	
	/**
	 * 获取row行col列元素的值
	 * @param row
	 * @param col
	 * @return
	 */
	double get(int row, int col);
	
	/**
	 * m(row,col) += value
	 * @param row
	 * @param col
	 * @param value
	 */
	void add(int row, int col,double value);
	
	/**
	 * get all non-zero element, instead of iterator
	 * 获取所有非零元素，不使用迭代子
	 * @return
	 */
	Map<Integer,Map<Integer,Double>> getAll();
	
	/**
	 * M(nRowBase + row, nColBase + col) = values in map
	 * 将参数map中的所有元素值赋值到本矩阵，其中map中所有的行号都加上nRowBase，
	 * 所有的列号都加上nColBase
	 * @param nRowBase
	 * @param nColBase
	 * @param map
	 */
	void setAll(int nRowBase, int nColBase, Map<Integer,Map<Integer,Double>> map);
	
	/**
	 * y = M*x
	 * 矩阵向量相乘，速度较慢，适用于快速代数方程组求解的矩阵接口参见AlgebraMatrix
	 * @param x
	 * @param y
	 */
	void mult(Vector x, Vector y);
	
	/**
	 * A=A'
	 */
	Matrix trans();
	
	/**
	 * Deep copy
	 * 深拷贝
	 * 
	 * @return
	 */
	Matrix copy();
	
	/**
	 * 清空矩阵中的所有元素
	 */
	void clear();
	
	/**
	 * 打印矩阵元素
	 */
	void print();
	
	Matrix setName(String name);
	String getName();
}