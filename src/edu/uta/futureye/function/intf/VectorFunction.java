package edu.uta.futureye.function.intf;

import java.util.List;
import java.util.Map;

import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.function.Variable;

public interface VectorFunction {
	/**
	 * Returns the value of vector function at <tt>v</tt>
	 * 返回向量值函数在<tt>v</tt>点的值
	 * 
	 * @param v
	 * @return
	 */
	Vector value(Variable v);
	
	/**
	 * Set variable names of the vector function
	 * 设置向量值函数的变量名
	 * 
	 * @param varNames
	 */
	void setVarNames(List<String> varNames);
	
	/**
	 * Returns variable names of the vector function
	 * 返回向量值函数的变量名列表
	 * 
	 * @return
	 */
	List<String> varNames();
	
	/**
	 * Get dimension of vector valued function
	 * 获得向量函数的维度
	 * 
	 * @return
	 */
	int getDim();
	
	/**
	 * Set dimension of vector valued function
	 * 设置向量函数的维度
	 * 
	 * @return
	 */
	void setDim(int dim);
	
	/**
	 * Set component <tt>index</tt> to function <tt>value</tt>
	 * 将分量<tt>index</tt>的值设置为<tt>value</tt>
	 * 
	 * @param index
	 * @param value
	 */
	void set(int index, Function value);
	
	/**
	 * Get function at <tt>index</tt>
	 * 获得分量<tt>index</tt>对应的函数
	 * 
	 * @param index
	 * @return
	 */
	Function get(int index);

	/**
	 * Composite vector valued function
	 * 向量值函数的符合函数
	 * 
	 * @param fInners
	 * @return
	 */
	VectorFunction compose(Map<String,Function> fInners);

	///////////////////////////////////////////////
	
	/**
	 * <code>fi(x)=gi(x), i=1...dim</code>
	 * 将向量函数<tt>\vec{g}(x)</tt>的值赋值给<tt>\vec{f}(x)</tt>
	 * 
	 * @param <code>\vec{g}(x)=(g1(x),g2(x),...,gn(x)</code>
	 */
	VectorFunction set(VectorFunction g);
	
	/**
	 * <code>fi(x)=a*gi(x), i=1...dim</code>
	 * 将向量函数<tt>a*g(x)</tt>的值赋值给<tt>f(x)</tt>
	 * 
	 * @param <code>\vec{g}(x)=(g1(x),g2(x),...,gn(x)</code>
	 */
	VectorFunction set(double a, VectorFunction g);
	
	/**
	 * <code>\vec{f}(x) = \vec{f}(x) + \vec{g}(x)</code>
	 * 
	 * @param <code>\vec{g}(x)=(g1(x), g2(x), ..., gn(x))</code>
	 */
	VectorFunction add(VectorFunction g);
	
	/**
	 * <code>\vec{f}(x) = \vec{f}(x) + a*\vec{g}(x)</code>
	 * 
	 * @param a
	 * @param <code>\vec{g}(x)=(g1(x), g2(x), ..., gn(x))</code>
	 */
	VectorFunction add(double a, VectorFunction g);
	
	/**
	 * <code>\vec{f}(x) = a*\vec{f}(x)</code>
	 * 
	 * @param a
	 * @return
	 */
	VectorFunction scale(double a);
	
	/**
	 * <code>\vec{f}(x) = a*\vec{f}(x)</code>
	 * 
	 * @param a
	 * @return
	 */
	VectorFunction ax(double a);
	
	/**
	 * <code>\vec{f}(x) = a*\vec{f}(x) + \vec{g}(x)</code>
	 * 
	 * @param a
	 * @param <code>\vec{g}(x)=(g1(x), g2(x), ..., gn(x))</code>
	 * @return
	 */
	VectorFunction axpy(double a, VectorFunction g);
	
	/**
	 * Dot product, returns 
	 * <code>f1(x)*g1(x) + f2(x)*g2(x) + ... + fn(x)*gn(x)</code>
	 * 点乘（内积）
	 * 
	 * @param <code>\vec{g}(x) = (g1(x), g2(x), ..., gn(x))</code>
	 * @return
	 */
	Function dot(VectorFunction g);
	
	/**
	 * Dot product, returns
	 * <code>f1(x)*g1 + f2(x)*g2 + ... + fn(x)*gn</code>
	 * 点乘（内积）
	 * 
	 * @param <code>\vec{g} = (g1, g2, ..., gn)</code>
	 * @return
	 */
	Function dot(Vector g);	
	
	////////////////////////////////////////////////////
	
	/**
	 *  Add
	 * (f1)   (g1)   (f1+g1)
	 * (f2) + (g2) = (f2+g2)
	 * (..)   (..)   ( ... )
	 * (fn)   (gn)   (fn+gn)
	 * 
	 * @param g
	 * @return
	 */
	VectorFunction A(VectorFunction g);
	VectorFunction A(Vector v);
	
	/**
	 *  Subtract
	 * (f1)   (g1)   (f1-g1)
	 * (f2) - (g2) = (f2-g2)
	 * (..)   (..)   ( ... )
	 * (fn)   (gn)   (fn-gn)
	 * 
	 * @param g
	 * @return
	 */
	VectorFunction S(VectorFunction g);
	VectorFunction S(Vector v);
	
	/**
	 *  Multiply (componentwise)
	 * (f1)   (g1)   (f1*g1)
	 * (f2) * (g2) = (f2*g2)
	 * (..)   (..)   ( ... )
	 * (fn)   (gn)   (fn*gn)
	 * 
	 * @param g
	 * @return
	 */
	VectorFunction M(VectorFunction g);
	VectorFunction M(Vector v);	
	
	/**
	 *  Divide (componentwise)
	 * (f1)   (g1)   (f1/g1)
	 * (f2) / (g2) = (f2/g2)
	 * (..)   (..)   ( ... )
	 * (fn)   (gn)   (fn/gn)
	 * 
	 * @param g
	 * @return
	 */
	VectorFunction D(VectorFunction g);
	VectorFunction D(Vector v);
	
	/////////////////////////////////////////////////
	
	/**
	 * Deep copy
	 * 深拷贝
	 * 
	 * @return
	 */
	VectorFunction copy();
	
	/**
	 * Print the component values of the vector function
	 */
	void print();
}
