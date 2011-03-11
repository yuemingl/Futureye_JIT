package edu.uta.futureye.function.intf;

import java.util.List;
import java.util.Map;

import edu.uta.futureye.function.Variable;

public interface Function {
	/**
	 * 返回函数值
	 * @param v
	 * @return
	 */
	double value(Variable v);
	
	/**
	 * 设置变量名称
	 * @param varNames
	 */
	void setVarNames(List<String> varNames);
	
	/**
	 * 返回所有变量名称
	 * @return
	 */
	List<String> varNames();
	
	/**
	 * this+f
	 * @param f
	 * @return
	 */
	Function P(Function f);
	
	/**
	 * this-f
	 * @param f
	 * @return
	 */
	Function M(Function f);
	
	/**
	 * this*f
	 * @param f
	 * @return
	 */
	Function X(Function f);
	
	/**
	 * this/f
	 * @param f
	 * @return
	 */
	Function D(Function f);
	
	/**
	 * 复合函数
	 * @param e.g. fInners: Map[ x = x(r,s), y = y(r,s) ]
	 * @return e.g. f = f(x,y) = f( x(r,s),y(r,s) )
	 */
	Function compose(Map<String,Function> fInners);
	
	/**
	 * f(x).d(x) := \frac{ \partial{this} }{ \partial{varName} }
	 * 
	 * 关于varName的一阶导数
	 * @param varName
	 * @return
	 */
	Function d(String varName);
	
	/**
	 * For constant function only
	 * @return
	 */
	double value();

}
