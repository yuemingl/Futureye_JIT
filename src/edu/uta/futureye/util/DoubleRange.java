package edu.uta.futureye.util;

/**
 * Range with double type of bound
 * 
 * @author liuyueming
 *
 */
public class DoubleRange {
	public double min;
	public double max;
	
	public DoubleRange(double min, double max) {
		this.min = min;
		this.max = max;
	}
	
	/**
	 * Return true if <tt>min < x < max</tt>
	 * @param x
	 * @return
	 */
	public boolean isInc(double x) {
		if(x>min && x<max) return true;
		else return false;
	}
	
	/**
	 * Return true if <tt>min <= x < max</tt>
	 * @param x
	 * @return
	 */
	public boolean isIncL(double x) {
		if(x>=min && x<max) return true;
		else return false;
	}	
	
	/**
	 * Return true if <tt>min < x <= max</tt>
	 * @param x
	 * @return
	 */
	public boolean isIncR(double x) {
		if(x>min && x<=max) return true;
		else return false;
	}
	
	/**
	 * Return true if <tt>min <= x <= max</tt>
	 * @param x
	 * @return
	 */
	public boolean isIncLR(double x) {
		if(x>=min && x<=max) return true;
		else return false;
	}
	
}
