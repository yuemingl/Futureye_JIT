/**
 * Copyright (c) 2010, nkliuyueming@gmail.com. All rights reserved.
 * 
 * 
 */
package edu.uta.futureye.function.intf;

import java.util.List;
import java.util.Map;

import com.sun.org.apache.bcel.internal.generic.ConstantPoolGen;
import com.sun.org.apache.bcel.internal.generic.InstructionFactory;
import com.sun.org.apache.bcel.internal.generic.InstructionHandle;
import com.sun.org.apache.bcel.internal.generic.InstructionList;
import com.sun.org.apache.bcel.internal.generic.MethodGen;

import edu.uta.futureye.bytecode.CompiledFunc;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.VariableArray;

/**
 * Mathematical function interface
 * 
 */
public interface MathFunc {
	/**
	 * Get the function name
	 * 
	 * @return
	 */
	String getName();
	
	/**
	 * Set a name for the function
	 * 
	 * @param name
	 * @return
	 */
	MathFunc setName(String name);
	
	/**
	 * Set free variable names of the function
	 * 
	 * @param varNames
	 * @return
	 */
	MathFunc setVarNames(List<String> varNames);
	
	/**
	 * Return all free variable names of the function
	 * <p>
	 * 返回所有自变量名称
	 * 
	 * @return
	 */
	List<String> getVarNames();
	
	/**
	 * Return function value at variable v
	 * <p>
	 * 返回自变量v对应的函数值
	 * 
	 * @param v
	 * @return double function value
	 */
	double apply(Variable v);
	
	double apply(double ...args);
	
	double apply(Element e, Node n, double ...args);
	
	/**
	 * Return function value at variable v with cache enabled
	 * <p>
	 * 返回自变量v对应的函数值，支持缓存。在一个复杂的函数表达式中，
	 * 对于同一函数多次求值的情况，使用缓存可以提高计算速度。例如在
	 * 有限元问题的数值积分中，Jacobian在被积函数中可能多处出现，
	 * 此时可以利用缓存对象在第一次求值后将其缓存起来，后续的求值调用
	 * 只需要从缓存中获取即可，不需要反复求值。
	 * 
	 * @param v
	 * @param cache
	 * @return
	 */
	double apply(Variable v, Map<Object,Object> cache);
	
	/**
	 * For more efficiency, this interface can be used to return function values
	 * at an array of variables at once.
	 * 
	 * @param valAry VariableArray object which represents an array of variable values
	 * @param cache Cache for efficient evaluation of functions. <tt>null</tt> parameter will disable the cache mechanism
	 * @return Function values evaluated at array of variables <tt>valAry</tt> 
	 */
	double[] applyAll(VariableArray valAry, Map<Object, Object> cache);
	
	/**
	 * Set the index of each variable in parameter 'args' for the 'apply' methods
	 *  
	 * @param argIdx
	 * @return
	 */
	MathFunc setArgIdx(Map<String, Integer> argsMap); 
	
	/**
	 * Get the index map of free variables of the function
	 * @return
	 */
	Map<String, Integer> getArgIdxMap();
	
	/**
	 * Add: this + g
	 */
	MathFunc A(MathFunc g);
	
	/**
	 * Add: this + g
	 */
	MathFunc A(double g);
	
	/**
	 * Subtract: this - g
	 */
	MathFunc S(MathFunc g);
	
	/**
	 * Subtract: this - g
	 */
	MathFunc S(double g);
	
	/**
	 * Multiply: this*g
	 */
	MathFunc M(MathFunc g);
	
	/**
	 * Multiply: this*g
	 */
	MathFunc M(double g);
	
	/**
	 * Divide: this/g
	 */
	MathFunc D(MathFunc g);
	
	/**
	 * Divide: this/g
	 */	
	MathFunc D(double g);
	
	/**
	 * Composition function 
	 * 复合函数
	 * <p><blockquote><pre>
	 * e.g.
	 * Function fx = FX.fx;
	 * Function fr = new FX("r");
	 * Function fOut = fx.M(fx).S(FC.c1);
	 * System.out.println(fOut); //x*x - 1.0
	 * 
	 * Function fIn = fr.M(fr).A(FC.c1);
	 * System.out.println(fIn); //r*r + 1.0
	 * 
	 * //Construct a map to define variable mapping
	 * Map<String,Function> map = new HashMap<String,Function>();
	 * map.put("x", fIn); //x=r*r + 1.0
	 * Function fComp = fOut.compose(map);
	 * System.out.println(fComp); //x(r)*x(r) - 1.0, where x=r*r + 1.0
	 * </pre></blockquote>
	 * 
	 * @param fInners: Variable map e.g.[ x = x(r,s), y = y(r,s) ]
	 * @return composed function e.g. f = f(x,y) = f( x(r,s),y(r,s) )
	 */
	MathFunc compose(Map<String,MathFunc> fInners);
	
	/**
	 * First derivative with respect to <code>varName</code> 
	 * <p>
	 * 关于<code>varName</code>的一阶导数：
	 * <p>
	 * f(x)._d(x) := \frac{ \partial{this} }{ \partial{varName} }
	 * 
	 * @param varName
	 * @return
	 */
	MathFunc diff(String varName);

	/**
	 * Returns true if it is a constant function
	 * 
	 * @return
	 */
	boolean isConstant();
	
	/**
	 * Shallow copy
	 * 
	 * @return
	 */
	MathFunc copy();
	
	/**
	* Return the expression (a string) of the function
	*
	* @return
	*/
	String getExpr();
	
	/**
	 * Definition of operand order
	 */
	static int OP_ORDER0 = 0; //Parenthesis first
	static int OP_ORDER1 = 1; //Exponents next
	static int OP_ORDER2 = 2; //Multiply and Divide next
	static int OP_ORDER3 = 3; //Add and Subtract last of all
	
	/**
	 * Get order of operations (Priority Rules for Arithmetic)
	 * <blockquote><pre>
	 * 0 Brackets first
	 * 1 Exponents next
	 * 2 Multiply and Divide next
	 * 3 Add and Subtract last of all.
	 * </blockquote></pre>  
	 */	
	int getOpOrder();
	
	/**
	 * Set order of operations (Priority Rules for Arithmetic)
	 * <blockquote><pre>
	 * 0 Brackets first
	 * 1 Exponents next
	 * 2 Multiply and Divide next
	 * 3 Add and Subtract last of all.
	 * </blockquote></pre>  
	 * @param order
	 */
	void setOpOrder(int order);
	
	/**
	 * Implement this function for your own compilation
	 * 
	 * @param clsName
	 * @param mg
	 * @param cp
	 * @param factory
	 * @param il
	 * @param argsMap Arguments name=>index pair
	 * @param argsStartPos Start index of 'args' in generated function apply()
	 * @param funcRefsMap An array of objects in the expression that implements MathFun 
	 * @return
	 */
	InstructionHandle bytecodeGen(String clsName, MethodGen mg, 
			ConstantPoolGen cp, InstructionFactory factory, 
			InstructionList il, Map<String, Integer> argsMap, 
			int argsStartPos, Map<MathFunc, Integer> funcRefsMap);

	/**
	 * Compile the function to bytecode with order of arguments from
	 * the call of getVarNames()
	 */
	CompiledFunc compile();
	
	/**
	 * Compile the function to bytecode with a specified order of
	 * arguments
	 * 
	 * @param varNames
	 * @return
	 */
	CompiledFunc compile(String[] varNames);
	
	//////////////Operator overloading support through Java-OO//////////////////
	/**
	 * Operator overloading support:
	 * 
	 * MathFunc a = 5;
	 * 
	 */
	MathFunc valueOf(int v);
	MathFunc valueOf(long v);
	MathFunc valueOf(float v) ;
	MathFunc valueOf(double v);
	
	/**
	 * Operator overload support:
	 * a+b
	 */
	MathFunc add(MathFunc other);
	MathFunc add(int other);
	MathFunc addRev(int other);
	MathFunc add(long other);
	MathFunc addRev(long other);
	MathFunc add(float other);
	MathFunc addRev(float other);
	MathFunc add(double other);
	MathFunc addRev(double other);
	
	/**
	 * Operator overload support:
	 * a-b
	 */
	MathFunc subtract(MathFunc other);
	MathFunc subtract(int other);
	MathFunc subtractRev(int other);
	MathFunc subtract(long other);
	MathFunc subtractRev(long other);
	MathFunc subtract(float other);
	MathFunc subtractRev(float other);
	MathFunc subtract(double other);
	MathFunc subtractRev(double other);
	
	/**
	 * Operator overload support:
	 * a*b
	 */
	MathFunc multiply(MathFunc other);
	MathFunc multiply(int other);
	MathFunc multiplyRev(int other);
	MathFunc multiply(long other);
	MathFunc multiplyRev(long other);
	MathFunc multiply(float other);
	MathFunc multiplyRev(float other);
	MathFunc multiply(double other);
	MathFunc multiplyRev(double other);
	
	/**
	 * Operator overload support:
	 * a/b
	 */
	MathFunc divide(MathFunc other);
	MathFunc divide(int other);
	MathFunc divideRev(int other);
	MathFunc divide(long other);
	MathFunc divideRev(long other);
	MathFunc divide(float other);
	MathFunc divideRev(float other);
	MathFunc divide(double other);
	MathFunc divideRev(double other);
	
	/**
	 * Operator overload support:
	 * -a
	 */
	MathFunc negate();

}
