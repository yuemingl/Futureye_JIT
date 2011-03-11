package edu.uta.futureye.algebra.intf;

/**
 * 为矩阵向量运算优化的向量接口
 * 
 * @author liuyueming
 *
 */
public interface AlgebraVector {
	/**
	 * 获得向量的维度
	 * @return
	 */
	public int getDim();
	
	/**
	 * 获得向量的数组
	 * @return
	 */
	public double[] getData();
	
	/**
	 * this = v
	 * @param v
	 */
	public AlgebraVector set(AlgebraVector v);

	/**
	 * this += v
	 * @param a
	 * @param v
	 */
	public AlgebraVector plus(AlgebraVector v);	

	/**
	 * this -= v
	 * @param a
	 * @param v
	 */
	public AlgebraVector minus(AlgebraVector v);	

	/**
	 * this += a*v
	 * @param a
	 * @param v
	 */
	public AlgebraVector add(double a, AlgebraVector v);	

	/**
	 * this *= a
	 * @param a
	 * @return
	 */
	public AlgebraVector scale(double a);

	/**
	 * this *= a
	 * 等价于 scale(double a);
	 * @param a
	 * @return
	 */
	public AlgebraVector ax(double a);
	
	/**
	 * x = a*x+y
	 * 等价于 add(double a, AlgebraVector v);
	 * @param a
	 * @param b
	 * @return
	 */
	public AlgebraVector axpy(double a, AlgebraVector y);
	
	/**
	 * x = (a*x).*y
	 * @param a
	 * @param y
	 * @return
	 */
	public AlgebraVector axmy(double a, AlgebraVector y);
	
	/**
	 * x.y
	 * 内积
	 * @param y
	 * @return
	 */
	public double dot(AlgebraVector y);
	
	/**
	 * 二范数
	 * @return
	 */
	public double norm2();
	
	/**
	 * 无穷范数
	 * @return
	 */
	public double normInf();
	
	/**
	 * print the component values of this vector
	 */
	public void print();

}
