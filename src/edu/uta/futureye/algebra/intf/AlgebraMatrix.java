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
	 * C = A*B
	 * 
	 * (A=this)
	 * 
	 * @param x
	 * @param y
	 */
	public void mult(AlgebraMatrix B, AlgebraMatrix C);
	
	/**
	 * Get A'
	 * @return A'
	 */
	public AlgebraMatrix getTrans();
	
	/**
	 * print the component values of this matrix
	 */
	public void print();

}
