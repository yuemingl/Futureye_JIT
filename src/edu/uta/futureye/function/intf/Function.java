package edu.uta.futureye.function.intf;

import java.util.List;
import java.util.Map;

import edu.uta.futureye.function.Variable;

/**
 * Function interface
 * 
 * @author liuyueming
 *
 */
public interface Function {
	/**
	 * Return function value at variable v
	 * 返回自变量v对应的函数值
	 * 
	 * @param v
	 * @return double function value
	 */
	double value(Variable v);
	
	/**
	 * Set function variable names
	 * 设置函数自变量名称，对于复合函数，只设置外层自变量名称
	 * 关于复合函数的构造 @see compose()
	 * 
	 * @param varNames
	 */
	void setVarNames(List<String> varNames);
	
	/**
	 * Return all variable names of the function
	 * 返回所有自变量名称
	 * 
	 * @return
	 */
	List<String> varNames();
	
	/**
	 * Add
	 * 
	 * @param g
	 * @return f+g, f==this
	 */
	Function A(Function g);
	Function A(double g);
	
	/**
	 * Subtract
	 * 
	 * @param g
	 * @return f-g, f==this
	 */
	Function S(Function g);
	Function S(double g);
	
	/**
	 * Multiply
	 * 
	 * @param g
	 * @return f*g, f==this
	 */
	Function M(Function g);
	Function M(double g);
	
	/**
	 * Divide
	 * 
	 * @param g
	 * @return f/g, f==this
	 */
	Function D(Function g);
	Function D(double g);
	
	/**
	 * Composition function
	 * 复合函数
	 * e.g.
	 *  Function fx = FX.fx;
	 *	Function fr = new FX("r");
	 *	Function fOut = fx.M(fx).S(FC.c1);
	 *	System.out.println(fOut); //x*x - 1.0
	 *	Function fIn = fr.M(fr).A(FC.c1);
	 *  System.out.println(fIn); //r*r + 1.0
	 *	//Construct a map to define variable mapping
	 *	Map<String,Function> map = new HashMap<String,Function>();
	 *	map.put("x", fIn); //x=r*r + 1.0
	 *	Function fComp = fOut.compose(map);
	 *	System.out.println(fComp); //x(r)*x(r) - 1.0, where x=r*r + 1.0
	 *
	 * @param e.g. fInners: Variable map[ x = x(r,s), y = y(r,s) ]
	 * @return e.g. f = f(x,y) = f( x(r,s),y(r,s) )
	 */
	Function compose(Map<String,Function> fInners);
	
	/**
	 * 关于varName的一阶导数：
	 * f(x)._d(x) := \frac{ \partial{this} }{ \partial{varName} }
	 * 
	 * @param varName
	 * @return
	 */
	Function _d(String varName);

	
	/**
	 * Return function (for constant function only)
	 * 
	 * @return
	 */
	double value();
	
	/**
	 * Deep copy
	 * @return
	 */
	Function copy();
	
	//////////////////For printing expression only///////////////////////////
	/**
	 * If function name is not null, the name is printed instead of the expression
	 */
	String getFName();
	void setFName(String name);
	
	/**
	 * Order of Operations (Priority Rules for Arithmetic)
	 *   0 Brackets first
	 *   1 Exponents next
	 *   2 Multiply and Divide next
	 *   3 Add and Subtract last of all.
	 */
	static int OP_ORDER0 = 0;
	static int OP_ORDER1 = 1;
	static int OP_ORDER2 = 2;
	static int OP_ORDER3 = 3;
	int getOpOrder();
	void setOpOrder(int order);
}
