package edu.uta.futureye.function.basic;

import java.util.HashMap;
import java.util.Map;

import edu.uta.futureye.function.AbstractFunction;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.intf.Function;

/**
 * Constant function: f=c
 * 
 * @author liuyueming
 *
 */
public class FC extends AbstractFunction{
	//Predefined constant
	public static FC c0 = new FC(0.0);
	public static FC c1 = new FC(1.0);
	public static FC PI = new FC(Math.PI);
	public static FC E = new FC(Math.E);
	
	protected double val = 0.0;
	protected static Map<Double, FC> cs = new HashMap<Double, FC>();
	
	public FC() {
	}
	
	public FC(double val) {
		this.val = val;
	}
	
	/**
	 * 返回常值函数，并且保存在静态Map中，以便多次使用，节省内存
	 * 注意：不要使用该函数生成大量常数，否则内存在程序退出前不会释放
	 * @param v
	 * @return
	 */
	public static FC c(double v) {
		FC c = cs.get(v);
		if(c == null) {
			c = new FC(v);
			cs.put(v, c);
			return c;
		} else {
			return c;
		}
	}
	
	@Override
	public double value(Variable v) {
		return val;
	}

	@Override
	public double value() {
		return val;
	}
	
	@Override
	public Function _d(String varName) {
		return c0;
	}
	
	@Override
	public int getOpOrder() {
		return OP_ORDER0;
	}
	
	@Override
	public Function copy() {
		return new FC(this.val);
	}
	
	@Override
	public String toString() {
		return String.valueOf(val);
	}
}
