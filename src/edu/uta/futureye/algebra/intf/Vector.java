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
	 * 设置向量的维度
	 * @param dim
	 */
	public void setDim(int dim);
	
	/**
	 * 获得向量的维度
	 * @return
	 */
	public int getDim();

	/**
	 * 将分量index的值设置为value
	 * @param index
	 * @param value
	 */
	public void set(int index, double value);
	
	/**
	 * 将向量v的值赋值给this
	 * @param v
	 */
	public void set(Vector v);

	/**
	 * 获得分量index的值
	 * @param index
	 * @return
	 */
	public double get(int index);

	/**
	 * v(index) += value
	 * 分量index的值再加上value
	 * @param index
	 * @param value
	 */
	public void add(int index,double value);
	
	/**
	 * this += a*v
	 * 
	 * @param a
	 * @param v
	 */
	public void add(double a, Vector v);

	/**
	 * 复制自己，深拷贝
	 * @return
	 */
	public Vector copy();
	
	/**
	 * 点乘（内积）
	 * Dot product
	 * 
	 * @param b = (b1 b2 ... bn)'
	 * @return = a1*b1 + a2*b2 + ... + an*bn
	 */
	public double dot(Vector b);
	
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
	 * 清空向量的所有元素
	 */
	public void clear();
	
	/**
	 * print the component values of this vector
	 * 打印向量元素
	 */
	public void print();
	

}
