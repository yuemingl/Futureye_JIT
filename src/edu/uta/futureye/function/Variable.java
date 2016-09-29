package edu.uta.futureye.function;

import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Map.Entry;

import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.geometry.Point;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.util.Constant;

/**
 * Function arguments (Independent variables of a function)
 * 
 * TODO 2011-4-26
 * 是否考虑定义Variable接口，这样可以实现多个类，
 * 1.只包括values
 * 2.只包括index
 * 3.values, index
 * 4.values, index, element
 * 这样的好处是可以减少Variable对象的大小
 * 
 * 
 * 
 * @author liuyueming
 *
 */
public class Variable {
	/**
	 * We use LinkedHashMap here. LinkedHashMap is hash table and linked list implementation of the Map interface, 
	 * with predictable iteration order.
	 * <p>
	 * 遍历时保证变量顺序
	*/ 
	protected Map<String,Double> valMap = new LinkedHashMap<String,Double>();

	//protected double[] valArray = new double[9];
	//protected boolean[] useArray = {false,false,false,false,false,false,false,false,false};
	//protected boolean bApplyRestirct = false;
	
	//Node Index
	protected int index = 0;
	
	//Element 变量中可以携带单元对象，见class DuDn
	protected Element element = null;
	
	public Variable() {
	}
	
	///////////////////////////////////////////////
	/**
	 * 构造一维自变量并赋值
	 */
	public Variable(double val) {
		//valArray[0] = val;
		//useArray[0] = true;
		valMap.put(Constant.x, val);
	}
	
	/**
	 * 获取一维自变量值
	 * @return
	 */
	public double get() {
		//return valArray[0];
		return valMap.values().iterator().next();
	}
	
	////////////////////////////////////////////////

	public Variable(String name, double val) {
//		int id = VN.getID(name);
//		if(id != -1) {
//			valArray[id] = val;
//			useArray[id] = true;
//		}
//		else
			valMap.put(name, val);
	}
	
	public Variable(VarPair first, VarPair ...pairs) {
//		int id = first.getVarID();
//		if(id != -1) {
//			valArray[id] = first.value;
//			useArray[id] = true;
//		} else
			valMap.put(first.name, first.value);
		for(int i=0;i<pairs.length;i++) {
//			id = pairs[i].getVarID();
//			if(id != -1) {
//				valArray[id] = pairs[i].value;
//				useArray[id] = true;
//			} else
				valMap.put(pairs[i].name, pairs[i].value);
		}
	}
	
	/**
	 * 返回自变量名名称对应的值
	 * @param name
	 * @return
	 */
	public double get(String name) {
//		int id = VN.getID(name);
//		if(id != -1)
//			return valArray[id];
//		else
			return valMap.get(name);
	}
	public double get(VN name) {
		//return valArray[name.getID()];
		return valMap.get(VN.names[name.getID()]);
	}
	
	/**
	 * Alias of get(String name), used in ScalaFEM as syntactic sugar: 
	 * <code>v(name)</code>
	 * 
	 * @param name
	 * @return <tt>v(name)</tt>
	 */
	public double apply(String name) {
		return valMap.get(name);
	}
	
	public double apply(VN name) {
		return valMap.get(name);
	}
	
	/**
	 * 
	 * @param pair
	 * @return
	 */
	public double get(VarPair pair) {
//		int id = pair.getVarID();
//		if(id != -1)
//			valArray[id] = pair.value;
//		else
			pair.value = valMap.get(pair.name);
		return pair.value;
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
//		int id = VN.getID(name);
//		if(id != -1) {
//			valArray[id] = val;
//			useArray[id] = true;
//		} else
			valMap.put(name, val);
		return this;
	}
	public Variable set(VN name, double val) {
//		int id = name.getID();
//		valArray[id] = val;
//		useArray[id] = true;
		valMap.put(VN.names[name.getID()],val);
		return this;
	}
	
	public Variable set(VarPair pair) {
//		int id = pair.getVarID();
//		if(id != -1) {
//			valArray[id] = pair.value;
//			useArray[id] = true;
//		} else
			valMap.put(pair.name, pair.value);
		return this;
	}

	public Map<String,Double> getNameValuePairs() {
//		for(int i=0;i<9;i++)
//			if(useArray[i]) valMap.put(VN.names[i], valArray[i]);
		return valMap;
	}

	public double[] getVarValues() {
		double[] vals = new double[this.valMap.size()];
		int i = 0;
		for(Entry<String, Double> e : valMap.entrySet()) {
			vals[i] = e.getValue();
			i++;
		}
		return vals;		
	}
	
	public String[] getVarNames() {
		String[] names = new String[this.valMap.size()];
		int i = 0;
		for(Entry<String, Double> e : valMap.entrySet()) {
			names[i] = e.getKey();
			i++;
		}
		return names;
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
	public static Variable createFrom(MathFunc fun, Point point, int index) {
		if(fun == null)
			return null;
		Variable var = new Variable();
		var.setIndex(index);
		if(fun.getVarNames() != null) {
			int ic = 1;
			for(String vn : fun.getVarNames()) {
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
		return this.getNameValuePairs().toString();
	}
	
	public static void main(String[] args) {
		Variable v1 = new Variable();
		System.out.println(v1);
		Variable v2 = new Variable(new VarPair("x",1.0));
		System.out.println(v2);
		Variable v3 = new Variable(new VarPair("x",1.0),
				new VarPair("y",2.0));
		System.out.println(v3);
		
		Variable v = new Variable("x",5.0).set("y",6.0);
		System.out.println(v.get("x"));
		System.out.println(v.get("y"));
		v.setIndex(20);
		System.out.println(v.getIndex());
		Element e = new Element();
		e.globalIndex = 100;
		v.setElement(e);
		System.out.println(v.getElement());
	}
	
}
