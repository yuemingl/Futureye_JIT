package edu.uta.futureye.function;

import java.util.LinkedHashMap;
import java.util.Map;

import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.geometry.Point;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.util.Constant;

/**
 * Function arguments (Independent variables of a function)
 * @author liuyueming
 *
 */
public class Variable {
	//LinkedHashMap 遍历时保证变量顺序
	protected Map<String,Double> values = new LinkedHashMap<String,Double>();

	//protected boolean bApplyRestirct = false;
	
	//Node Index
	protected int index = 0;
	
	//变量中可以携带element，见class DuDn
	protected Element element = null;
	
	public Variable() {
	}
	
	///////////////////////////////////////////////
	/**
	 * 构造一维自变量并赋值
	 */
	public Variable(double val) {
		values.put(Constant.x, val);
	}
	
	/**
	 * 获取一维自变量值
	 * @return
	 */
	public double get() {
		//TODO Need check values.size==1 ?
		return values.values().iterator().next();
	}
	
	////////////////////////////////////////////////

	public Variable(String name, double val) {
		values.put(name, val);
	}
	
	public Variable(VarPair fitst, VarPair ...pairs) {
		values.put(fitst.name, fitst.value);
		for(int i=0;i<pairs.length;i++)
			values.put(pairs[i].name, pairs[i].value);
	}
	
	/**
	 * 返回自变量名名称对应的值
	 * @param name
	 * @return
	 */
	public double get(String name) {
		return values.get(name);
	}
	
	/**
	 * 设置自变量的值：名称、值对， 返回变量自身，以便链式表达：
	 * 例如：设置二维自变量x=1,y=2：
	 * Variable v = new Variable("x",1).set("y",2);
	 * 
	 * @param name
	 * @param val
	 * @return
	 */
	public Variable set(String name, double val) {
		values.put(name, val);
		return this;
	}

	public Map<String,Double> getValues() {
		return values;
	}

	
//	public void applyRestirct(boolean flag) {
//		bApplyRestirct = flag;
//	}
//	
//	/**
//	 * 如果函数对自变量有限制，例如二维形函数限制到一维边界上，isRestrict()会提示函数求值时是否应用该限制
//	 * @return
//	 */
//	public boolean isRestrict() {
//		return bApplyRestirct;
//	}

	
	/**
	 * 根据fun的自变量个数和Point的值（以及index，如果需要的话），
	 * 创建一个Variable的对象
	 * @param fun
	 * @param point
	 * @param index
	 * @return
	 */
	public static Variable createFrom(Function fun, Point point, int index) {
		if(fun == null)
			return null;
		Variable var = new Variable(index);
		if(fun.varNames() != null) {
			int ic = 1;
			for(String vn : fun.varNames()) {
				var.set(vn, point.coord(ic));
				ic++;
			}
		} else {
			//VectorBasedFunction
		}
		return var;
	}
	
	/**
	 * 变量中携带向量索引（编号）
	 * @param index
	 */
	public Variable setIndex(int index) {
		this.index = index;
		return this;
	}
	public int getIndex() {
		return index;
	}
	
	/**
	 * 变量中携带单元对象
	 * @return
	 */
	public Element getElement() {
		return element;
	}
	public Variable setElement(Element element) {
		this.element = element;
		return this;
	}
	
	public String toString() {
		return values.toString();
	}
	
	public static void main(String[] args) {
		Variable v1 = new Variable();
		System.out.println(v1);
		Variable v2 = new Variable(new VarPair("x",1.0));
		System.out.println(v2);
		Variable v3 = new Variable(new VarPair("x",1.0),
				new VarPair("y",2.0));
		System.out.println(v3);
	}
	
}
