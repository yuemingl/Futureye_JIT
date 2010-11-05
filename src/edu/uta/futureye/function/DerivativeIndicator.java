package edu.uta.futureye.function;

import java.util.LinkedHashMap;
import java.util.Map;

/**
 * 导数指示器，指明关于那个变量的几阶导数
 * @author liuyueming
 *
 */
public class DerivativeIndicator {
	protected Map<String,Integer> ind = new LinkedHashMap<String,Integer>();
	
	public DerivativeIndicator() {
	}
	
	/**
	 * e.g. "x_2,y_1" = fun_xx_y
	 * @param arg
	 * @return
	 */
	public DerivativeIndicator(String arg) {
		//TODO
		for(int i=0;i<arg.length();i+=2) {
			set(String.valueOf(arg.charAt(i)),
					Integer.parseInt(String.valueOf(arg.charAt(i+1))));
		}
	}
	
	public void set(String varName, int degree) {
		ind.put(varName, degree);
	}
	
	public Map<String,Integer> get() {
		return ind;
	}
	
	public String toString() {
		return ind.toString();
	}
}
