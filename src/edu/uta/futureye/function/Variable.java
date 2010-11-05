package edu.uta.futureye.function;

import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Map.Entry;

import edu.uta.futureye.core.intf.Point;
import edu.uta.futureye.function.intf.Function;

/**
 * Function arguments (Independent variables of a function)
 * @author liuyueming
 *
 */
public class Variable {
	protected Map<String,Double> values = new LinkedHashMap<String,Double>();
//	protected boolean bApplyRestirct = false;
	//Node Index
	protected int index = 0;
	
	public Variable() {
	}
	public Variable(int index) {
		this.index = index;
	}	
	public Variable(String name, double val) {
		values.put(name, val);
	}
	
	public Variable(VarPair fitst, VarPair ...pairs) {
		values.put(fitst.name, fitst.value);
		for(int i=0;i<pairs.length;i++)
			values.put(pairs[i].name, pairs[i].value);
	}
	
	public double get(String name) {
		return values.get(name);
	}
	
	public Variable(double val) {
		values.put("x", val);
	}
	/**
	 * 适用于一维自变量
	 */
	public double get() {
		return values.values().iterator().next();
	}	
	
	public void set(String name, double val) {
		values.put(name, val);
	}

	public Map<String,Double> getValues() {
		return values;
	}

	public String toString() {
		return values.toString();
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
	
	
	public void setIndex(int index) {
		this.index = index;
	}
	
	public int getIndex() {
		return index;
	}
	
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
