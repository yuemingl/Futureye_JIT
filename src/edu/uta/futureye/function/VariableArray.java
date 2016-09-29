package edu.uta.futureye.function;

import java.util.LinkedHashMap;
import java.util.Map;

import edu.uta.futureye.core.Element;
import edu.uta.futureye.util.Constant;

/**
 * Array of Function arguments. 
 * In mathematical context, this class represents an array of 
 * independent variables of functions.
 * 
 * @author liuyueming
 */
public class VariableArray {
	private int length = 0;
	//LinkedHashMap 遍历时保证变量顺序
	protected Map<String,double[]> valMap = new LinkedHashMap<String,double[]>();

	//Node Indices
	protected int[] indices = null;
	
	//Element 变量中可以携带单元对象，见class DuDn
	protected Element element = null;
	
	/**
	 * Construct an empty variable array. 
	 * Use set() to set values for a variable of function
	 */
	public VariableArray() {
	}
	
	/**
	 * Construct values array of 1D variable x
	 * <p>
	 * 构造一维自变量x的值数组
	 */
	public VariableArray(double[] valAry) {
		length = valAry.length;
		valMap.put(Constant.x, valAry);
	}
	
	/**
	 * Construct values array of 1D variable <tt>name</tt>
	 * 
	 * @param name Variable name
	 * @param valAry Array of variable values
	 */
	public VariableArray(String name, double[] valAry) {
		length = valAry.length;
		valMap.put(name, valAry);
	}
	
	/**
	 * Get values array of 1D variable x
	 * <p>
	 * 获取一维自变量x的值数组
	 * @return
	 */
	public double[] get() {
		return valMap.get(Constant.x);
	}
	
	/**
	 * 返回自变量名名称对应的值
	 * @param name
	 * @return
	 */
	public double[] get(String name) {
		return valMap.get(name);
	}
	
	/**
	 * Set value array for variable <tt>name</tt>
	 * 
	 * @param name variable name
	 * @param valAry value array 
	 * @return this variable
	 */
	public VariableArray set(String name, double[] valAry) {
		length = valAry.length;
		valMap.put(name, valAry);
		return this;
	}
	
	/**
	 * Get the map of variable name and value array
	 * 
	 * @return
	 */
	public Map<String,double[]> getValues() {
		return valMap;
	}

	/**
	 * Set an integer index for this variable. 
	 * This index can be an index of a vector.
	 * <p>
	 * 变量中携带向量索引（编号）数组
	 * 
	 * @param index
	 * @return this variable
	 */
	public VariableArray setIndices(int[] index) {
		this.indices = index;
		return this;
	}
	
	/**
	 * Get the integer index which is set by <tt>setIndex()</tt>
	 * 
	 * @return
	 */
	public int[] getIndices() {
		return indices;
	}
	

	/**
	 * Set an element reference for this variable.
	 * <p>
	 * 设置变量中携带的单元对象
	 * @param element finite element object
	 * @return this variable
	 */
	public VariableArray setElement(Element element) {
		this.element = element;
		return this;
	}
	
	/**
	 * Get the element object which is set by <tt>setElement()</tt>
	 * <p>
	 * 获得变量中携带的单元对象
	 * @return
	 */
	public Element getElement() {
		return element;
	}
	
	public int length() {
		return this.length;
	}
	
	public String toString() {
		return this.getValues().toString();
	}
	
}
