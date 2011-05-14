package edu.uta.futureye.algebra.intf;

/**
 * General vector interface
 * 一般向量接口，不适用于求解代数方程组（求解代数方程组的向量接口参见AlgebraVector）
 * 用途：
 * 1. 二维、三维“空间向量”，用于向量函数
 * 2. 用来保存右端向量合成步骤中的局部和全局向量
 * 3. 矩阵稀疏存贮的“解向量”，用于打印输出
 * 
 * @author liuyueming
 * 
 */
public interface Vector {
	/**
	 * Get dimension of the vector
	 * 获得向量的维度
	 * 
	 * @return
	 */
	int getDim();
	
	/**
	 * Set dimension of the vector
	 * 设置向量的维度
	 * 
	 * @param dim
	 */
	void setDim(int dim);
	
	/**
	 * Get <code>x(index)</code>
	 * 获得分量<tt>index</tt>的值
	 * 
	 * @param index
	 * @return <code>x(index)</code>
	 */
	double get(int index);

	/**
	 * <code>x(index) = value</code>
	 * 将分量<tt>index</tt>的值设置为<tt>value</tt>
	 * 
	 * @param index
	 * @param value
	 */
	void set(int index, double value);
	
	/**
	 * <code>x(index) += value</code>
	 * 分量<tt>index</tt>的值再加上<tt>value</tt>，
	 * 适用于装配刚度矩阵和载荷向量
	 * 
	 * @param index
	 * @param value
	 */
	void add(int index,double value);
	
	////////////////////////////////////////////////
	
	/**
	 * <code>x(i)=y(i), i=1...dim</code>
	 * 将向量y的值赋值给x
	 * 
	 * @param y
	 */
	Vector set(Vector y);
	
	/**
	 * <code>x(i)=a*y(i), i=1...dim</code>
	 * 将向量a*y的值赋值给x
	 * 
	 * @param y
	 */
    Vector set(double a, Vector y);
	
	/**
	 * <code>x = x + y</code>
	 * 
	 * @param a
	 * @param y
	 */
	Vector add(Vector y);
	
	/**
	 * <code>x = x + a*y</code>
	 * 
	 * @param a
	 * @param y
	 */
	Vector add(double a, Vector y);
	
	/**
	 * <code>x = a*x</code>
	 * 
	 * @param a
	 * @return
	 */
	Vector scale(double a);
	
	/**
	 * <code>x = a*x</code>
	 * 
	 * @param a
	 * @return
	 */
	Vector ax(double a);
	
	/**
	 * <code>x = a*x + y</code>
	 * 
	 * @param a
	 * @param y
	 * @return
	 */
	Vector axpy(double a, Vector y);
	
	/**
	 * xi = a*xi*yi
	 * @param a
	 * @param y
	 * @return
	 */
	Vector axMuly(double a, Vector y);
	
	/**
	 * xi = a*xi/yi
	 * @param a
	 * @param y
	 * @return
	 */
	Vector axDivy(double a, Vector y);
	
	/**
	 * <code>x = x + dv</code>
	 * 
	 * @param dv
	 * @return
	 */
	Vector shift(double dv);
	
	/**
	 * Dot product, returns <code>x1*y1 + x2*y2 + ... + xn*yn</code>
	 * 点乘（内积）
	 * 
	 * @param y = (y1 y2 ... yn)'
	 * @return
	 */
	double dot(Vector y);
	
	/**
	 * <code>sum(abs(x))</code>
	 * 一范数
	 * 
	 * @return
	 */
	double norm1();

	/**
	 * <code>sqrt(sum(abs(x).^2))</code>
	 * 二范数
	 * 
	 * @return
	 */
	double norm2();
	
	/**
	 * <code>max(abs(x))</code>
	 * 无穷范数
	 * 
	 * @return
	 */
	double normInf();
	
	///////////////////////////////////////////////////////
	
	/**
	 * Deep copy
	 * 深拷贝
	 * 
	 * @return
	 */
	Vector copy();

	/**
	 * Set vector dimension to zero. (SparseVector only?) 
	 * 清空向量的所有元素（是否定义稀疏向量接口？）
	 */
	void clear();
	
	/**
	 * Print the component values of the vector
	 * 打印向量元素
	 */
	void print();
}
