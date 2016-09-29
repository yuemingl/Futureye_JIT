package edu.uta.futureye.algebra.intf;



/**
 * <blockquote><pre>
 * General matrix interface
 * 一般矩阵接口，不适用于求解代数方程组（求解代数方程组的矩阵接口参见AlgebraMatrix）
 * 用途：
 * 1. 用来保存刚度矩阵合成步骤中的局部和全局矩阵
 * 2. 打印输出矩阵
 * </blockquote></pre>
 * 
 * @author liuyueming
 *
 */
public interface Matrix {
	/**
	 * 矩阵零元素阈值
	 */
	public static double zeroEps = 1e-15;
	
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
	 * 获取row行col列元素的值
	 * @param row
	 * @param col
	 * @return <tt>m(row,col)</tt>
	 */
	double get(int row, int col);
	
	/**
	 * Alias of get(int row, int col), used in ScalaFEM as syntactic sugar: 
	 * <code>val a = m(row,col)</code>
	 * 
	 * @param row
	 * @param col
	 * @return <tt>m(row,col)</tt>
	 */
	double apply(int row, int col);

	/**
	 * 设置row行col列元素的值
	 * @param row
	 * @param col
	 * @param value
	 */
	void set(int row, int col, double value);
	
	/**
	 * Alias of set(int row, int col, double value), used in ScalaFEM as syntactic sugar: 
	 * <code>m(row,col)=value</code>
	 * 
	 * @param row
	 * @param col
	 * @param value
	 */
	void update(int row, int col, double value);

	/**
	 * m(row,col) += value
	 * @param row
	 * @param col
	 * @param value
	 */
	void add(int row, int col,double value);
	
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
	 * Print matrix information
	 * <p>
	 * 打印矩阵元素
	 */
	void print();
	
	/**
	 * Set matrix name for printing purpose or using in Matlab as variable name
	 * 
	 * @param name Matrix name
	 * @return <tt>this</tt> for convenience only
	 */
	Matrix setName(String name);
	
	/**
	 * Get matrix name
	 * 
	 * @return Matrix name
	 */
	String getName();
	
	/**
	 * Write this matrix to file <tt>fileName</tt> with Matlab mat file format
	 * <p>
	 * If more than one matrix need to be written in a single mat file use <tt>MatlabMatFileWriter</tt> instead.
	 * @param fileName
	 */
	void writeMatFile(String fileName);
	
	/**
	 * Write this matrix to file <tt>fileName</tt> with simple text file format
	 * @param fileName
	 */
	void writeSimpleFile(String fileName);
	
	//writeMatrixMarketFile(String fileName);
}