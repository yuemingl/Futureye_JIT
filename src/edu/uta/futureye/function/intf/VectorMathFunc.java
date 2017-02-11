package edu.uta.futureye.function.intf;

import java.util.List;
import java.util.Map;

import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.VariableArray;

public interface VectorMathFunc {
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
	VectorMathFunc compose(Map<String,MathFunc> fInners);

	///////////////////////////////////////////////
	
	/**
	 * <code>fi(x)=gi(x), i=1...dim</code>
	 * <p>
	 * 将向量函数<tt>\vec{g}(x)</tt>的值赋值给<tt>\vec{f}(x)</tt>
	 * 
	 * @param <code>\vec{g}(x)=(g1(x),g2(x),...,gn(x)</code>
	 */
	VectorMathFunc set(VectorMathFunc g);
	
	/**
	 * <code>fi(x)=a*gi(x), i=1...dim</code>
	 * <p>
	 * 将向量函数<tt>a*g(x)</tt>的值赋值给<tt>f(x)</tt>
	 * 
	 * @param <code>\vec{g}(x)=(g1(x),g2(x),...,gn(x)</code>
	 */
	VectorMathFunc set(double a, VectorMathFunc g);
	
	/**
	 * <code>\vec{f}(x) = \vec{f}(x) + \vec{g}(x)</code>
	 * 
	 * @param <code>\vec{g}(x)=(g1(x), g2(x), ..., gn(x))</code>
	 */
	//VectorMathFunc add(VectorMathFunc g);
	
	/**
	 * <code>\vec{f}(x) = \vec{f}(x) + a*\vec{g}(x)</code>
	 * 
	 * @param a
	 * @param <code>\vec{g}(x)=(g1(x), g2(x), ..., gn(x))</code>
	 */
	VectorMathFunc add(double a, VectorMathFunc g);
	
	/**
	 * <code>\vec{f}(x) = a*\vec{f}(x)</code>
	 * 
	 * @param a
	 * @return
	 */
	VectorMathFunc scale(double a);
	
	/**
	 * <code>\vec{f}(x) = a*\vec{f}(x)</code>
	 * 
	 * @param a
	 * @return
	 */
	VectorMathFunc ax(double a);
	
	/**
	 * <code>\vec{f}(x) = a*\vec{f}(x) + \vec{g}(x)</code>
	 * 
	 * @param a
	 * @param <code>\vec{g}(x)=(g1(x), g2(x), ..., gn(x))</code>
	 * @return
	 */
	VectorMathFunc axpy(double a, VectorMathFunc g);
	
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
	MathFunc dot(VectorMathFunc g);
	
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
	VectorMathFunc A(VectorMathFunc g);
	
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
	VectorMathFunc A(Vector v);
	
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
	VectorMathFunc S(VectorMathFunc g);
	
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
	VectorMathFunc S(Vector v);
	
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
	VectorMathFunc M(VectorMathFunc g);
	
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
	VectorMathFunc M(Vector v);	
	
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
	VectorMathFunc D(VectorMathFunc g);
	
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
	VectorMathFunc D(Vector v);
	
	/////////////////////////////////////////////////
	
	/**
	 * Deep copy
	 * 深拷贝
	 * 
	 * @return
	 */
	VectorMathFunc copy();
	
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
	VectorMathFunc setFName(String name);
	
	//////////////Operator overloading support through Java-OO//////////////////
	/**
	 * Operator overloading support:
	 * 
	 * MathFunc a = 5;
	 * 
	 */
	VectorMathFunc valueOf(int v);
	VectorMathFunc valueOf(long v);
	VectorMathFunc valueOf(float v) ;
	VectorMathFunc valueOf(double v);
	
	/**
	 * Operator overload support:
	 * a+b
	 */
	VectorMathFunc add(VectorMathFunc other);
	VectorMathFunc add(int other);
	VectorMathFunc addRev(int other);
	VectorMathFunc add(long other);
	VectorMathFunc addRev(long other);
	VectorMathFunc add(float other);
	VectorMathFunc addRev(float other);
	VectorMathFunc add(double other);
	VectorMathFunc addRev(double other);
	
	/**
	 * Operator overload support:
	 * a-b
	 */
	VectorMathFunc subtract(VectorMathFunc other);
	VectorMathFunc subtract(int other);
	VectorMathFunc subtractRev(int other);
	VectorMathFunc subtract(long other);
	VectorMathFunc subtractRev(long other);
	VectorMathFunc subtract(float other);
	VectorMathFunc subtractRev(float other);
	VectorMathFunc subtract(double other);
	VectorMathFunc subtractRev(double other);
	
	/**
	 * Operator overload support:
	 * a*b
	 */
	VectorMathFunc multiply(VectorMathFunc other);
	VectorMathFunc multiply(int other);
	VectorMathFunc multiplyRev(int other);
	VectorMathFunc multiply(long other);
	VectorMathFunc multiplyRev(long other);
	VectorMathFunc multiply(float other);
	VectorMathFunc multiplyRev(float other);
	VectorMathFunc multiply(double other);
	VectorMathFunc multiplyRev(double other);
	
	/**
	 * Operator overload support:
	 * a/b
	 */
	VectorMathFunc divide(VectorMathFunc other);
	VectorMathFunc divide(int other);
	VectorMathFunc divideRev(int other);
	VectorMathFunc divide(long other);
	VectorMathFunc divideRev(long other);
	VectorMathFunc divide(float other);
	VectorMathFunc divideRev(float other);
	VectorMathFunc divide(double other);
	VectorMathFunc divideRev(double other);
	
	/**
	 * Operator overload support:
	 * -a
	 */
	VectorMathFunc negate();
	
}
