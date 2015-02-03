package edu.uta.futureye.algebra.intf;

/**
 * <blockquote><pre>
 * Vector interface designed for fast algebra operation
 * 为矩阵向量运算优化的向量接口
 * 
 * </blockquote></pre>
 * 
 * @author liuyueming
 *
 */
public interface AlgebraVector {
	/**
	 * Get the dimension of the vector
	 * 获得向量的维度
	 * 
	 * @return
	 */
	public int getDim();
	
	/**
	 * Get <tt>double</tt> array of the vector 
	 * 获得向量的数组
	 * 
	 * @return
	 */
	public double[] getData();
	
	/**
	 * <code>x = y</code>
	 * 
	 * @param y
	 */
	public AlgebraVector set(AlgebraVector y);

	/**
	 * <code>x = a*y</code>
	 * 
	 * @param y
	 */
	public AlgebraVector set(double a, AlgebraVector y);
	
	/**
	 * <code>x = x + y</code>
	 * 
	 * @param a
	 * @param y
	 */
	public AlgebraVector add(AlgebraVector y);
	
	/**
	 * <code>x = x - y</code>
	 * 
	 * @param a
	 * @param y
	 */
	public AlgebraVector subtract(AlgebraVector y);	

	/**
	 * <code>x = x + a*y</code>
	 * 
	 * @param a
	 * @param y
	 */
	public AlgebraVector add(double a, AlgebraVector y);

	/**
	 * <code>x = a*x</code>
	 * 
	 * @param a
	 * @return
	 */
	public AlgebraVector scale(double a);

	/**
	 * <code>x = a*x</code>
	 * Alias of <code>scale(double a)</code>
	 * <p>
	 * 等价于 <code>scale(double a)</code>
	 * 
	 * @param a
	 * @return
	 */
	public AlgebraVector ax(double a);
	
	/**
	 * <code>x = a*x+y</code>
	 * <p>
	 * Notice: Different from <code>add(double a, AlgebraVector y)</code>
	 * <p>
	 * 注意：与<code>add(double a, AlgebraVector y)</code>有区别
	 * 
	 * @param a
	 * @param b
	 * @return
	 */
	public AlgebraVector axpy(double a, AlgebraVector y);
	
	/**
	 * <code>x = (a*x).*y</code>
	 * 
	 * @param a
	 * @param y
	 * @return
	 */
	public AlgebraVector axmy(double a, AlgebraVector y);
	
	/**
	 * Inner product <code>x.y</code>
	 * <p>
	 * 内积
	 * 
	 * @param y
	 * @return
	 */
	public double dot(AlgebraVector y);
	
	/**
	 * 1 Norm
	 * 
	 * @return
	 */
	public double norm1();
	
	/**
	 * 2 Norm
	 * 
	 * @return
	 */
	public double norm2();
	
	/**
	 * Infinity norm (Maximum norm)
	 * 
	 * @return
	 */
	public double normInf();
	
	/**
	 * Print vector entries
	 */
	public void print();

}
