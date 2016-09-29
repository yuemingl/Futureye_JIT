package edu.uta.futureye.function.intf;

import java.util.List;
import java.util.Map;

import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.VariableArray;

public interface VectorFunction {
	/**
	 * Returns the value of vector function at <tt>v</tt>
	 * <p>
	 * 返回向量值函数在<tt>v</tt>点的值
	 * 
	 * @param v
	 * @return
	 */
	Vector value(Variable v);
	
	/**
	 * Returns the array of value of vector function at <tt>valAry</tt>
	 * 
	 * @param v
	 * @param cache
	 * @return
	 */
	Vector[] valueArray(VariableArray valAry, Map<Object,Object> cache);
	
	/**
	 * Set variable names of the vector function
	 * <p>
	 * 设置向量值函数的变量名
	 * 
	 * @param varNames
	 */
	void setVarNames(List<String> varNames);
	
	/**
	 * Returns variable names of the vector function
	 * <p>
	 * 返回向量值函数的变量名列表
	 * 
	 * @return
	 */
	List<String> varNames();
	
	/**
	 * Get dimension of vector valued function
	 * <p>
	 * 获得向量函数的维度
	 * 
	 * @return
	 */
	int getDim();
	
	/**
	 * Set dimension of vector valued function
	 * <p>
	 * 设置向量函数的维度
	 * 
	 * @return
	 */
	void setDim(int dim);
	
	/**
	 * Set component <tt>index</tt> to function <tt>value</tt>
	 * <p>
	 * 将分量<tt>index</tt>的值设置为<tt>value</tt>
	 * 
	 * @param index
	 * @param value
	 */
	void set(int index, MathFunc value);
	
	/**
	 * Get function at <tt>index</tt>
	 * <p>
	 * 获得分量<tt>index</tt>对应的函数
	 * 
	 * @param index
	 * @return
	 */
	MathFunc get(int index);

	/**
	 * Composite vector valued function
	 * <p>
	 * 向量值函数的符合函数
	 * 
	 * @param fInners
	 * @return
	 */
	VectorFunction compose(Map<String,MathFunc> fInners);

	///////////////////////////////////////////////
	
	/**
	 * <code>fi(x)=gi(x), i=1...dim</code>
	 * <p>
	 * 将向量函数<tt>\vec{g}(x)</tt>的值赋值给<tt>\vec{f}(x)</tt>
	 * 
	 * @param <code>\vec{g}(x)=(g1(x),g2(x),...,gn(x)</code>
	 */
	VectorFunction set(VectorFunction g);
	
	/**
	 * <code>fi(x)=a*gi(x), i=1...dim</code>
	 * <p>
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
	 * <p>
	 * <code>f1(x)*g1(x) + f2(x)*g2(x) + ... + fn(x)*gn(x)</code>
	 * <p>
	 * 点乘（内积）
	 * 
	 * @param <code>\vec{g}(x) = (g1(x), g2(x), ..., gn(x))</code>
	 * @return
	 */
	MathFunc dot(VectorFunction g);
	
	/**
	 * Dot product, returns
	 * <p>
	 * <code>f1(x)*g1 + f2(x)*g2 + ... + fn(x)*gn</code>
	 * <p>
	 * 点乘（内积）
	 * 
	 * @param <code>\vec{g} = (g1, g2, ..., gn)</code>
	 * @return
	 */
	MathFunc dot(Vector g);	
	
	////////////////////////////////////////////////////
	
	/**
	 * Add
	 * <p><blockquote><pre>
	 * (f1)   (g1)   (f1+g1)
	 * (f2) + (g2) = (f2+g2)
	 * (..)   (..)   ( ... )
	 * (fn)   (gn)   (fn+gn)
	 * </blockquote></pre>
	 * 
	 * @param g
	 * @return
	 */
	VectorFunction A(VectorFunction g);
	
	/**
	 * Add
	 * <p><blockquote><pre>
	 * (f1)   (v1)   (f1+v1)
	 * (f2) + (v2) = (f2+v2)
	 * (..)   (..)   ( ... )
	 * (fn)   (vn)   (fn+vn)
	 * </blockquote></pre>
	 * 
	 * @param v
	 * @return
	 */
	VectorFunction A(Vector v);
	
	/**
	 * Subtract
	 * <p><blockquote><pre>
	 * (f1)   (g1)   (f1-g1)
	 * (f2) - (g2) = (f2-g2)
	 * (..)   (..)   ( ... )
	 * (fn)   (gn)   (fn-gn)
	 * </blockquote></pre>
	 * 
	 * @param g
	 * @return
	 */
	VectorFunction S(VectorFunction g);
	
	/**
	 *  Subtract
	 *  <p><blockquote><pre>
	 * (f1)   (v1)   (f1-v1)
	 * (f2) - (v2) = (f2-v2)
	 * (..)   (..)   ( ... )
	 * (fn)   (vn)   (fn-vn)
	 * </blockquote></pre>
	 * 
	 * @param v
	 * @return
	 */
	VectorFunction S(Vector v);
	
	/**
	 *  Multiply (componentwise) with vector function
	 *  <p><blockquote><pre>
	 * (f1)   (g1)   (f1*g1)
	 * (f2) * (g2) = (f2*g2)
	 * (..)   (..)   ( ... )
	 * (fn)   (gn)   (fn*gn)
	 * </blockquote></pre>
	 * 
	 * @param g
	 * @return
	 */
	VectorFunction M(VectorFunction g);
	
	/**
	 * Multiply (componentwise) with vector
	 *  <p><blockquote><pre>
	 * (f1)   (v1)   (v1*f1)
	 * (f2) * (v2) = (v2*f2)
	 * (..)   (..)   ( ... )
	 * (fn)   (vn)   (vn*fn)
	 * </blockquote></pre>
	 * 
	 * @param v
	 * @return
	 */
	VectorFunction M(Vector v);	
	
	/**
	 *  Divide (componentwise) by vector function
	 * <p><blockquote><pre>
	 * (f1)   (g1)   (f1/g1)
	 * (f2) / (g2) = (f2/g2)
	 * (..)   (..)   ( ... )
	 * (fn)   (gn)   (fn/gn)
	 * </blockquote></pre>
	 * 
	 * @param g
	 * @return
	 */
	VectorFunction D(VectorFunction g);
	
	/**
	 * Divide (componentwise) by vector
	 * <p><blockquote><pre>
	 * (f1)   (v1)   (f1/v1)
	 * (f2) / (v2) = (f2/v2)
	 * (..)   (..)   ( ... )
	 * (fn)   (vn)   (fn/vn)
	 * </blockquote></pre>
	 * 
	 * @param v: a Vector
	 * @return
	 */
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
	 * return the expression of function
	 * 
	 * @return
	 */
	String getExpression();
	
	/**
	 * Get function name
	 * <p>
	 * If function name is not null, the name instead of the expression of
	 * function is returned by <code>toString()</code> method
	 */
	String getFName();
	
	/**
	 * Set function name
	 * <p>
	 * If function name is not null, the name instead of the expression of
	 * function is returned by <code>toString()</code> method
	 * @param name
	 */
	VectorFunction setFName(String name);
	
}
