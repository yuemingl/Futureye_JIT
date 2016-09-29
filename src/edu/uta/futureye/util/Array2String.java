package edu.uta.futureye.util;

public class Array2String {
	public static String convert(int[] ary) {
		StringBuilder sb = new StringBuilder();
		sb.append("[");
		for(int i=0;i<ary.length-1;i++) {
			sb.append(ary[i]).append(",");
		}
		sb.append(ary[ary.length-1]).append("]");
		return sb.toString();
	}
	
	public static String convert(Integer[] ary) {
		StringBuilder sb = new StringBuilder();
		sb.append("[");
		for(int i=0;i<ary.length-1;i++) {
			sb.append(ary[i]).append(",");
		}
		sb.append(ary[ary.length-1]).append("]");
		return sb.toString();
	}
	
	public static String convert(double[] ary) {
		StringBuilder sb = new StringBuilder();
		sb.append("[");
		for(int i=0;i<ary.length-1;i++) {
			sb.append(ary[i]).append(",");
		}
		sb.append(ary[ary.length-1]).append("]");
		return sb.toString();
	}
	
	public static String convert(Double[] ary) {
		StringBuilder sb = new StringBuilder();
		sb.append("[");
		for(int i=0;i<ary.length-1;i++) {
			sb.append(ary[i]).append(",");
		}
		sb.append(ary[ary.length-1]).append("]");
		return sb.toString();
	}
}
