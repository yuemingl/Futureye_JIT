package edu.uta.futureye.application;

import edu.uta.futureye.function.AbstractMathFunc;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.util.FutureyeException;
import static edu.uta.futureye.function.operator.FMath.*;

public class ModelParam {
	static double bk = 0.1; //Background of mu_a = 0.1

	/**
	 * Get background value of mu_a (0.1)
	 * 
	 * @return
	 */
	public static MathFunc getBkMu_a() {
		return C(bk);
	}
	/**
	 * Generate mu_a for testing
	 * Background of mu_a = 0.1
	 * 
	 * type=1: one inclusion
	 * type=2: two inclusion
	 * @param incX
	 * @param incY
	 * @param incR
	 * @param maxMu_a
	 * @param type
	 */
	public static MathFunc getMu_a(double incX, double incY, double incR, 
			double maxMu_a, int type) {
		if(type == 1)
			return getMu_a(incX, incY, incR, maxMu_a, 1, 0.0, false);
		else
			return getMu_a(incX, incY, incR, maxMu_a, 2, 0.8, false);
	}
	
	/**
	 * Advanced function for creating mu_a
	 * 
	 * @param incX
	 * @param incY
	 * @param incR
	 * @param maxMu_a
	 * @param num_inclusion
	 * @param distance
	 * @param uniform
	 * @return
	 */
	public static MathFunc getMu_a(double incX, double incY, double incR, 
			double maxMu_a,int num_inclusion, double distance, boolean uniform) {
		final double fcx = incX;
		final double fcy = incY;
		final double fcr = incR;
		final double fmu_a = maxMu_a;
		final double dis = distance;
		final boolean funi = uniform;
		
		
		MathFunc rltMu_a = null;
		if(num_inclusion == 1) {
			rltMu_a = new AbstractMathFunc("x","y"){
				@Override
				public double apply(Variable v) {
					double dx = v.get("x")-fcx;
					double dy = v.get("y")-fcy;
					if(Math.sqrt(dx*dx+dy*dy) < fcr) {
						double r = fmu_a;
						if(!funi) {
							double theta = Math.PI/2 - (Math.PI/2)*Math.sqrt(dx*dx+dy*dy)/fcr;
							r = fmu_a*( Math.sin(Math.sin(Math.sin(
										Math.sin(Math.sin(Math.sin(
										Math.sin(Math.sin(Math.sin(theta
												)))))))))/0.502171); 
						}
						return r<bk?bk:r;
					}
					else
						return bk;
				}
			};
		} else if(num_inclusion == 2){
			rltMu_a = new AbstractMathFunc("x","y"){
				@Override
				public double apply(Variable v) {
					double dx = v.get("x")-fcx;
					double dy = v.get("y")-fcy;
					double dx1 = v.get("x")-(fcx+dis);
					double dy1 = v.get("y")-fcy;
					if(Math.sqrt(dx*dx+dy*dy) < fcr) {
						double r = fmu_a;
						if(!funi)
							r = fmu_a*Math.cos((Math.PI/2)*Math.sqrt(dx*dx+dy*dy)/fcr); 
						return r<bk?bk:r;
					}
					else if(Math.sqrt(dx1*dx1+dy1*dy1) < fcr) {
						double r = fmu_a;
						if(!funi)
							r = fmu_a*Math.cos((Math.PI/2)*Math.sqrt(dx1*dx1+dy1*dy1)/fcr); 
						return r<bk?bk:r;
					}
					else
						return bk;
				}
			};	
		} else if(num_inclusion == 3) { //Three inclusions
			rltMu_a = new AbstractMathFunc("x","y"){
				@Override
				public double apply(Variable v) {
					double dx = v.get("x")-fcx;
					double dy = v.get("y")-fcy;
					double dx1 = v.get("x")-(fcx+dis*Math.cos(Math.PI/3));
					double dy1 = v.get("y")-(fcy-dis*Math.sin(Math.PI/3));
					double dx2 = v.get("x")-(fcx-dis*Math.cos(Math.PI/3));
					double dy2 = v.get("y")-(fcy-dis*Math.sin(Math.PI/3));
					if(Math.sqrt(dx*dx+dy*dy) < fcr) {
						double r = fmu_a;
						if(!funi)
							r = fmu_a*Math.cos((Math.PI/2)*Math.sqrt(dx*dx+dy*dy)/fcr); 
						return r<bk?bk:r;
					}
					else if(Math.sqrt(dx1*dx1+dy1*dy1) < fcr) {
						double r = fmu_a;
						if(!funi)
							r = fmu_a*Math.cos((Math.PI/2)*Math.sqrt(dx1*dx1+dy1*dy1)/fcr); 
						return r<bk?bk:r;
					}
					else if(Math.sqrt(dx2*dx2+dy2*dy2) < fcr) {
						double r = fmu_a;
						if(!funi)
							r = fmu_a*Math.cos((Math.PI/2)*Math.sqrt(dx2*dx2+dy2*dy2)/fcr); 
						return r<bk?bk:r;
					}
					else
						return bk;
				}
			};			
		} else if(num_inclusion == 5) { //Five inclusions
			rltMu_a = new AbstractMathFunc("x","y"){
				@Override
				public double apply(Variable v) {
					double dx = v.get("x")-fcx;
					double dy = v.get("y")-fcy;
					double a = 36*Math.PI/180;
					double b = 72*Math.PI/180;
					double dx1 = v.get("x")-(fcx+dis*Math.cos(a));
					double dy1 = v.get("y")-(fcy-dis*Math.sin(a));
					double dx2 = v.get("x")-(fcx-dis*Math.cos(a));
					double dy2 = v.get("y")-(fcy-dis*Math.sin(a));
					double dx3 = v.get("x")-(fcx+dis*Math.cos(a)-dis*Math.cos(b));
					double dy3 = v.get("y")-(fcy-dis*Math.sin(a)-dis*Math.sin(b));
					double dx4 = v.get("x")-(fcx-dis*Math.cos(a)+dis*Math.cos(b));
					double dy4 = v.get("y")-(fcy-dis*Math.sin(a)-dis*Math.sin(b));
					if(Math.sqrt(dx*dx+dy*dy) < fcr) {
						double r = fmu_a;
						if(!funi)
							r = fmu_a*Math.cos((Math.PI/2)*Math.sqrt(dx*dx+dy*dy)/fcr); 
						return r<bk?bk:r;
					}
					else if(Math.sqrt(dx1*dx1+dy1*dy1) < fcr) {
						double r = fmu_a;
						if(!funi)
							r = fmu_a*Math.cos((Math.PI/2)*Math.sqrt(dx1*dx1+dy1*dy1)/fcr); 
						return r<bk?bk:r;
					}
					else if(Math.sqrt(dx2*dx2+dy2*dy2) < fcr) {
						double r = fmu_a;
						if(!funi)
							r = fmu_a*Math.cos((Math.PI/2)*Math.sqrt(dx2*dx2+dy2*dy2)/fcr); 
						return r<bk?bk:r;
					}
					else if(Math.sqrt(dx3*dx3+dy3*dy3) < fcr) {
						double r = fmu_a;
						if(!funi)
							r = fmu_a*Math.cos((Math.PI/2)*Math.sqrt(dx3*dx3+dy3*dy3)/fcr); 
						return r<bk?bk:r;
					}
					else if(Math.sqrt(dx4*dx4+dy4*dy4) < fcr) {
						double r = fmu_a;
						if(!funi)
							r = fmu_a*Math.cos((Math.PI/2)*Math.sqrt(dx4*dx4+dy4*dy4)/fcr); 
						return r<bk?bk:r;
					}
					else
						return bk;
				}
			};					
		} else {
			throw new FutureyeException("Unsupported parameter: num_inclusion="+num_inclusion);
		}
		return rltMu_a;
	}	
}
