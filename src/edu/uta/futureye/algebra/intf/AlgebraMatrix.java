package edu.uta.futureye.algebra.intf;

/**
 * 为矩阵向量运算优化的矩阵接口
 * 
 * @author liuyueming
 *
 */
public interface AlgebraMatrix {
	public int getRowDim();
	public int getColDim();
	
	/**
	 * y=A*x
	 * @param x
	 * @param y
	 */
	public void mult(AlgebraVector x, AlgebraVector y);

	/**
	 * print the component values of this matrix
	 */
	public void print();

}
