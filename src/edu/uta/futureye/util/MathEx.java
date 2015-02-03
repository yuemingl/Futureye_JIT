package edu.uta.futureye.util;

public class MathEx {
	public static double coth(double z) {
		double ez = Math.exp(z);
		double e_z = Math.exp(-z);
		return (ez+ez)/(ez-e_z);
	}
}
