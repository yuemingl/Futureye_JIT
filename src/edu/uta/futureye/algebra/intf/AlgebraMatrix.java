package edu.uta.futureye.algebra.intf;

/**
 * <blockquote><pre>
 * Matrix interface designed for fast algebra operation
 * 为矩阵向量运算优化的矩阵接口
 * 
 * </blockquote></pre>
 * 
 * @author liuyueming
 *
 */
public interface AlgebraMatrix {
	/**
	 * Get number of rows
	 */
	public int getRowDim();
	
	/**
	 * Get number of columns
	 * @return
	 */
	public int getColDim();
	
	/**
	 * Matrix vector multiplication
	 * y=this*x (y=A*x)
	 * 
	 * @param x
	 * @param y
	 */
	public void mult(AlgebraVector x, AlgebraVector y);

	/**
	 * Matrix matrix multiplication
	 * C = this*B (C=A*B)
	 * 
	 * @param x
	 * @param y
	 */
	public void mult(AlgebraMatrix B, AlgebraMatrix C);
	
	/**
	 * Get transpose of A
	 * 
	 * @return A'
	 */
	public AlgebraMatrix getTrans();
	
	/**
	 * print matrix entries
	 */
	public void print();

}
