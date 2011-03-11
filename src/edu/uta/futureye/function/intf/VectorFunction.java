package edu.uta.futureye.function.intf;

import java.util.List;

import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.function.Variable;

public interface VectorFunction {
	Vector value(Variable v);
	
	void setVarNames(List<String> varNames);
	List<String> varNames();
	
	/**
	 * 获得向量的维度
	 * @return
	 */
	int getDim();
	
	/**
	 * 将分量index的值设置为value
	 * @param index
	 * @param value
	 */
	void set(int index, Function value);
	
	/**
	 * 获得分量index的值
	 * @param index
	 * @return
	 */
	Function get(int index);

	/**
	 * 复制自己
	 * @return
	 */
	VectorFunction copy();
	
	/**
	 * 二范数
	 * @return
	 */
	Function norm2();
	
	/**
	 * 无穷范数
	 * @return
	 */
	Function normInf();
	
	/**
	 * 点乘（内积）
	 * Dot product
	 * 
	 * @param b = (b1(x) b2(x) ... bn(x))'
	 * @return = a1(x)*b1(x) + a2(x)*b2(x) + ... + an(x)*bn(x)
	 */
	Function dot(VectorFunction b);
	
	/**
	 * 点乘（内积）
	 * Dot product
	 * 
	 * @param b = (b1 b2 ... bn)'
	 * @return = a1(x)*b1 + a2(x)*b2 + ... + an(x)*bn
	 */
	Function dot(Vector b);
	
	/**
	 * print the component values of this vector function
	 */
	void print();
}
