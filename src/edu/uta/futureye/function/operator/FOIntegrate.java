package edu.uta.futureye.function.operator;

import java.util.HashMap;

import edu.uta.futureye.bytecode.CompiledFunc;
import edu.uta.futureye.function.VN;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.VariableArray;
import edu.uta.futureye.function.intf.MathFunc;

/**
 * Function Operator FOIntegrate: 
 *    Integrate functions on reference element
 * 
 * @author liuyueming
 *
 */
public class FOIntegrate{
	
	static double[] a2 = {
		-0.577350269189626,
		 0.577350269189626
		};
	static double[] h2 = {
		 1.0,
		 1.0
	};
	
	static double[] a3 = {
		-0.774596669241483,
		 0.0,
		 0.774596669241483
		 };
	static double[] h3 = {
		 0.555555555555556,
		 0.888888888888889,
		 0.555555555555556
	};
	
	static double[] a4 = {
		-0.861136311594953,
		-0.339981043584856,
		 0.339981043584856,
		 0.861136311594953
		 };
	static double[] h4 = {
		 0.347854845137454,
		 0.652145154862546,
		 0.652145154862546,
		 0.347854845137454
	};
	
	//Dean Joseph E. Flaherty (fea6.pdf)
	//Christoffel weights(h5) and roots(a5)
	static double[] a5 = {
		 0,
		 0.538469310105683,
		-0.538469310105683,
		 0.906179845938664,
		-0.906179845938664
	};
	static double[] h5 = {
		 0.568888888888889,
		 0.478628670499366,
		 0.478628670499366,
		 0.236926885056189,
		 0.236926885056189
		};	
	
	
    static double[] triW = {
        	0.06296959,
	        0.06619708,
	        0.06296959,
	        0.06619708,
	        0.06296959,
	        0.06619708,
	        0.11250000
        };
	static double[] triR = {
	        0.10128651,
	        0.47014206,
	        0.79742699,
	        0.47014206,
	        0.10128651,
	        0.05971587,
	        0.33333333
        };
	static double[] triS = {
			0.10128651,
			0.05971587,
			0.10128651,
			0.47014206,
			0.79742699,
			0.47014206,
			0.33333333
	    };
	
	/**
	 * 在参考三角形单元上积分
	 * @param integrand
	 * @param order
	 * @return Function
	 */	
	@SuppressWarnings("deprecation")
	public static double intOnTriangleRefElement(MathFunc integrand, int order) {
		double rlt = 0.0;
		if(order == 2) {
			Variable v = new Variable();
			v.set(VN.r, 0.333333333333333);
			v.set(VN.s, 0.333333333333333);
			v.set(VN.t, 0.333333333333333);
			rlt = 0.5*integrand.apply(v);
		} else if(order == 3) {
			
			Variable v1 = new Variable();
			Variable v2 = new Variable();
			Variable v3 = new Variable();
			
			v1.set(VN.r, 0.5); v1.set(VN.s, 0.5); v1.set(VN.t, 0.0);
			v2.set(VN.r, 0.0); v2.set(VN.s, 0.5); v2.set(VN.t, 0.5);
			v3.set(VN.r, 0.5); v3.set(VN.s, 0.0); v3.set(VN.t, 0.5);
			double pv1 = integrand.apply(v1);
			double pv2 = integrand.apply(v2);
			double pv3 = integrand.apply(v3);
			rlt = 0.5*0.333333333333333 * (pv1+pv2+pv3);
		} else if(order == 4) {
			double w123 = 25.0/48.0;
			double w4 = -27.0/48.0;
			
			Variable v1 = new Variable();
			Variable v2 = new Variable();
			Variable v3 = new Variable();
			Variable v4 = new Variable();
			
			v1.set(VN.r, 0.6); v1.set(VN.s, 0.2); v1.set(VN.t, 0.2);
			v2.set(VN.r, 0.2); v2.set(VN.s, 0.6); v2.set(VN.t, 0.2);
			v3.set(VN.r, 0.2); v3.set(VN.s, 0.2); v3.set(VN.t, 0.6);
			v4.set(VN.r, 0.333333333333333); v4.set(VN.s, 0.333333333333333); v4.set(VN.t, 0.333333333333333);
			
			rlt = w123*(integrand.apply(v1)+integrand.apply(v2)+integrand.apply(v3)) + w4*integrand.apply(v4);
			rlt /= 2;
		} else if(order == 5) {
			Variable v = new Variable();
			for(int i=0;i<7;i++) {
				v.set(VN.r, triR[i]); 
				v.set(VN.s, triS[i]); 
				v.set(VN.t, 1.0-triR[i]-triS[i]);
				rlt += triW[i]*integrand.apply(v);
			}
		} else {
			System.out.println("ERROR: intOnTriangleRefElement() Not supported order = "+order);
		}
		return rlt;
	}
	
	
	/**
	 * Integrate on 1D line segment reference element [-1,1]
	 * 
	 * @param integrand
	 * @param order
	 * @return
	 */
	@SuppressWarnings("deprecation")
	public static double intOnLinearRefElement(MathFunc integrand, int order) {
		double a1_1 = 0.0;
		double h1_1 = 2.0;
		
		double a2_1 = - 0.577350269189626;
		double a2_2 = - a2_1;
		double h2_1 = 1.0;
		double h2_2 = h2_1;
		
		double a3_1 = - 0.774596669241483;
		double a3_2 = 0.0;
		double a3_3 = - a3_1;
		double h3_1 = 0.555555555555556;
		double h3_2 = 0.888888888888889;
		double h3_3 = h3_1;
		
		double a4_1 = - 0.861136311594953;
		double a4_2 = - 0.339981043584856;
		double a4_3 = - a4_2;
		double a4_4 = - a4_1;
		double h4_1 = 0.347854845137454;
		double h4_2 = 0.652145154862546;
		double h4_3 = h4_2;
		double h4_4 = h4_1;
		
		Variable v = new Variable();
		double rlt = 0.0;
		if(order == 1) {
			v.set(VN.r, a1_1);
			rlt += h1_1*integrand.apply(v);
		} else if(order == 2) {
			v.set(VN.r, a2_1);
			rlt += h2_1*integrand.apply(v);
			v.set(VN.r, a2_2);
			rlt += h2_2*integrand.apply(v);
		} else if(order == 3) {
			v.set(VN.r, a3_1);
			rlt += h3_1*integrand.apply(v);
			v.set(VN.r, a3_2);
			rlt += h3_2*integrand.apply(v);
			v.set(VN.r, a3_3);
			rlt += h3_3*integrand.apply(v);
		} else if(order == 4) {
			v.set(VN.r, a4_1);
			rlt += h4_1*integrand.apply(v);
			v.set(VN.r, a4_2);
			rlt += h4_2*integrand.apply(v);
			v.set(VN.r, a4_3);
			rlt += h4_3*integrand.apply(v);
			v.set(VN.r, a4_4);
			rlt += h4_4*integrand.apply(v);
		} else if(order == 5) {
			for(int i=0;i<order;i++) {
				v.set(VN.r, a5[i]);
				rlt += h5[i]*integrand.apply(v);
			}
		} else {
			System.out.println("ERROR: intOnLinearRefElement() Not supported order = "+order);
		}
		return rlt;
	}
	
	/**
	 * Integrate on 2D rectangle reference element [-1,1]*[-1,1]
	 * 
	 * @param integrand
	 * @param order
	 * @return
	 */
	@SuppressWarnings("deprecation")
	public static double intOnRectangleRefElement(MathFunc integrand, int order) {
		double a1_1 = 0.0;
		double h1_1 = 4.0;
		double a2 = 0.577350269189626;
		//double h2 = 1.0;
		
		Variable v = new Variable();
		double rlt = 0.0;
		if(order == 1) {
			v.set(VN.r, a1_1);
			v.set(VN.s, a1_1);
			rlt += h1_1*integrand.apply(v);
		} else if(order == 2) {
			v.set(VN.r, a2);
			v.set(VN.s, a2);
			rlt += integrand.apply(v);
			v.set(VN.r, -a2);
			v.set(VN.s, a2);
			rlt += integrand.apply(v);
			v.set(VN.r, a2);
			v.set(VN.s, -a2);
			rlt += integrand.apply(v);
			v.set(VN.r, -a2);
			v.set(VN.s, -a2);
			rlt += integrand.apply(v);
		} else if(order == 5) {
			
//			for(int i=0;i<order;i++) {
//				for(int j=0;j<order;j++) {
//					v.set(VN.r, a5[i]);
//					v.set(VN.s, a5[j]);
//					rlt += h5[i]*h5[j]*integrand.value(v);
//				}
//			}
			
			VariableArray valAry = new VariableArray();
			double []ra = new double[25];
			double []sa = new double[25];
			double []wa = new double[25];
			int c = 0;
			for(int i=0;i<order;i++) {
				for(int j=0;j<order;j++) {
					ra[c] = a5[i];
					sa[c] = a5[j];
					wa[c] = h5[i]*h5[j];
					c++;
				}
			}
			valAry.set("r", ra);
			valAry.set("s", sa);
			double[] rltAry = integrand.applyAll(valAry,new HashMap<Object, Object>());
			//bugfix 2/3/2012 乘以wa[i]
			for(int i=0;i<rltAry.length;i++) 
				rlt += wa[i]*rltAry[i];
		} else {
			System.out.println("ERROR: intOnRectangleRefElement() Not supported order = "+order);
		}
		
		return rlt;
	}		
	
	/**
	 * 
	 * @param integrand
	 * @param order
	 * @return
	 */
	@SuppressWarnings("deprecation")
	public static double intOnTetrahedraRefElement(MathFunc integrand, int order) {
		double a1_1 = 0.25;
		double h1_1 = 1;
		
		double []a2 = {0.585410196624969,0.138196601125011,0.138196601125011,0.138196601125011};
		double h2 = 0.25;
		int [][]M24 = {{0,1,2,3},{1,0,2,3},{1,2,0,3},{1,2,3,0}};
		
		Variable v = new Variable();
		double rlt = 0.0;
		if(order == 1) {
			v.set(VN.r, a1_1);
			v.set(VN.s, a1_1);
			v.set(VN.t, a1_1);
			v.set(VN.u, a1_1);
			rlt += h1_1*integrand.apply(v);
		} else if(order ==2) {
			for(int i=0;i<M24.length;i++) {
				v.set(VN.r, a2[M24[i][0]]);
				v.set(VN.s, a2[M24[i][1]]);
				v.set(VN.t, a2[M24[i][2]]);
				v.set("u", a2[M24[i][3]]);
				rlt += h2*integrand.apply(v);
			}
		} else {
			System.out.println("ERROR: intOnTetrahedraRefElement() Not supported order = "+order);
		}
		
		return rlt;
	}
	
	public static double intOnHexahedraRefElement(MathFunc integrand, int order) {
//		Variable v = new Variable();
		VariableArray valAry = new VariableArray();
		double rlt = 0.0;
		if(order == 2) {
//			for(int i=0;i<order;i++) {
//				for(int j=0;j<order;j++) {
//					for(int k=0;k<order;k++) {
//						v.set(VN.r, a2[i]);
//						v.set(VN.s, a2[j]);
//						v.set(VN.t, a2[k]);
//						rlt += h2[i]*h2[j]*h2[k]*integrand.value(v);
//					}
//				}
//			}
			double []ra = new double[8];
			double []sa = new double[8];
			double []ta = new double[8];
			double []wa = new double[8];
			int c = 0;
			for(int i=0;i<order;i++) {
				for(int j=0;j<order;j++) {
					for(int k=0;k<order;k++) {
						ra[c] = a2[i];
						sa[c] = a2[j];
						ta[c] = a2[k];
						wa[c] = h2[i]*h2[j]*h2[k];
						c++;
					}
				}
			}
			valAry.set("r", ra);
			valAry.set("s", sa);
			valAry.set("t", ta);
			double[] rltAry = integrand.applyAll(valAry,new HashMap<Object, Object>());
			for(int i=0;i<rltAry.length;i++) 
				rlt += wa[i]*rltAry[i];

		} else if(order == 3) {
//			for(int i=0;i<order;i++) {
//				for(int j=0;j<order;j++) {
//					for(int k=0;k<order;k++) {
//						v.set(VN.r, a3[i]);
//						v.set(VN.s, a3[j]);
//						v.set(VN.t, a3[k]);
//						rlt += h3[i]*h3[j]*h3[k]*integrand.value(v);
//					}
//				}
//			}
			double []ra = new double[27];
			double []sa = new double[27];
			double []ta = new double[27];
			double []wa = new double[27];
			int c = 0;
			for(int i=0;i<order;i++) {
				for(int j=0;j<order;j++) {
					for(int k=0;k<order;k++) {
						ra[c] = a3[i];
						sa[c] = a3[j];
						ta[c] = a3[k];
						wa[c] = h3[i]*h3[j]*h3[k];
						c++;
					}
				}
			}
			valAry.set("r", ra);
			valAry.set("s", sa);
			valAry.set("t", ta);
			double[] rltAry = integrand.applyAll(valAry,new HashMap<Object, Object>());
			for(int i=0;i<rltAry.length;i++) 
				rlt += wa[i]*rltAry[i];
			
		} else if(order == 4) {
//			for(int i=0;i<order;i++) {
//				for(int j=0;j<order;j++) {
//					for(int k=0;k<order;k++) {
//						v.set(VN.r, a4[i]);
//						v.set(VN.s, a4[j]);
//						v.set(VN.t, a4[k]);
//						rlt += h4[i]*h4[j]*h4[k]*integrand.value(v);
//					}
//				}
//			}
			double []ra = new double[64];
			double []sa = new double[64];
			double []ta = new double[64];
			double []wa = new double[64];
			int c = 0;
			for(int i=0;i<order;i++) {
				for(int j=0;j<order;j++) {
					for(int k=0;k<order;k++) {
						ra[c] = a4[i];
						sa[c] = a4[j];
						ta[c] = a4[k];
						wa[c] = h4[i]*h4[j]*h4[k];
						c++;
					}
				}
			}
			valAry.set("r", ra);
			valAry.set("s", sa);
			valAry.set("t", ta);
			double[] rltAry = integrand.applyAll(valAry,new HashMap<Object, Object>());
			for(int i=0;i<rltAry.length;i++) 
				rlt += wa[i]*rltAry[i];
		} else if(order == 5) {
//			for(int i=0;i<order;i++) {
//				for(int j=0;j<order;j++) {
//					for(int k=0;k<order;k++) {
//						v.set(VN.r, a5[i]);
//						v.set(VN.s, a5[j]);
//						v.set(VN.t, a5[k]);
//						rlt += h5[i]*h5[j]*h5[k]*integrand.value(v);
//					}
//				}
//			}
			double []ra = new double[125];
			double []sa = new double[125];
			double []ta = new double[125];
			double []wa = new double[125];
			int c = 0;
			for(int i=0;i<order;i++) {
				for(int j=0;j<order;j++) {
					for(int k=0;k<order;k++) {
						ra[c] = a5[i];
						sa[c] = a5[j];
						ta[c] = a5[k];
						wa[c] = h5[i]*h5[j]*h5[k];
						c++;
					}
				}
			}
			valAry.set("r", ra);
			valAry.set("s", sa);
			valAry.set("t", ta);
			double[] rltAry = integrand.applyAll(valAry,new HashMap<Object, Object>());
			for(int i=0;i<rltAry.length;i++) 
				rlt += wa[i]*rltAry[i];
		}
		return rlt;
	}
	
	public static double intOnTriangleRefElement(CompiledFunc integrand, double[] params, int paramsStart, int order) {
		double rlt = 0.0;
		if(order == 2) {
			params[paramsStart] = 0.333333333333333;
			params[paramsStart+1] = 0.333333333333333;
			params[paramsStart+2] = 0.333333333333333;
			rlt = 0.5*integrand.apply(params);
		} else if(order == 3) {
			params[paramsStart] = 0.5; params[paramsStart+1] = 0.5; params[paramsStart+2] = 0.0; 
			double pv1 = integrand.apply(params);
			params[paramsStart] = 0.0; params[paramsStart+1] = 0.5; params[paramsStart+2] = 0.5; 
			double pv2 = integrand.apply(params);
			params[paramsStart] = 0.5; params[paramsStart+1] = 0.0; params[paramsStart+2] = 0.5; 
			double pv3 = integrand.apply(params);
			rlt = 0.5*0.333333333333333*(pv1 + pv2 + pv3);
		} else if(order == 4) {
			double w123 = 25.0/48.0;
			double w4 = -27.0/48.0;
			
			params[paramsStart] = 0.6; params[paramsStart+1] = 0.2; params[paramsStart+2] = 0.2; 
			double pv1 = integrand.apply(params);
			params[paramsStart] = 0.2; params[paramsStart+1] = 0.6; params[paramsStart+2] = 0.2; 
			double pv2 = integrand.apply(params);
			params[paramsStart] = 0.2; params[paramsStart+1] = 0.2; params[paramsStart+2] = 0.6; 
			double pv3 = integrand.apply(params);
			params[paramsStart] = 0.333333333333333; params[paramsStart+1] = 0.333333333333333; params[paramsStart+2] = 0.333333333333333; 
			double pv4 = 0.5*integrand.apply(params);
			
			rlt = 0.5*w123*(pv1 + pv2 + pv3) + w4*pv4;
		} else if(order == 5) {
			for(int i=0;i<7;i++) {
				params[paramsStart]   = triR[i]; 
				params[paramsStart+1] = triS[i]; 
				params[paramsStart+2] = 1.0-triR[i]-triS[i];
				rlt += triW[i]*integrand.apply(params);
			}
		}
		return rlt;
	}
	
	/**
	 * Integrate on 2D rectangle reference element [-1,1]*[-1,1]
	 */
	public static double intOnRectangleRefElement(CompiledFunc integrand, double[] params, int paramsStart, int order) {
		double a1_1 = 0.0;
		double h1_1 = 4.0;
		double a2 = 0.577350269189626;
		//double h2 = 1.0;
		
		double rlt = 0.0;
		if(order == 1) {
			params[paramsStart] = a1_1;
			params[paramsStart+1] = a1_1;
			rlt += h1_1*integrand.apply(params);
		} else if(order == 2) {
			params[paramsStart] = a2; 
			params[paramsStart+1] = a2;
			rlt += integrand.apply(params);
			params[paramsStart] = -a2; 
			params[paramsStart+1] = a2;
			rlt += integrand.apply(params);
			params[paramsStart] = a2; 
			params[paramsStart+1] = -a2;
			rlt += integrand.apply(params);
			params[paramsStart] = -a2; 
			params[paramsStart+1] = -a2;
			rlt += integrand.apply(params);
		} else if(order == 5) {
			for(int i=0;i<order;i++) {
				for(int j=0;j<order;j++) {
					params[paramsStart] = a5[i]; 
					params[paramsStart+1] = a5[j];
					rlt += h5[i]*h5[j]*integrand.apply(params);
				}
			}
		} else {
			System.out.println("ERROR: intOnRectangleRefElement() Not supported order = "+order);
		}
		
		return rlt;
	}

	public static double intOnLinearRefElement(CompiledFunc integrand, double[] params, int paramsStart, int order) {
		double a1_1 = 0.0;
		double h1_1 = 2.0;
		
		double a2_1 = - 0.577350269189626;
		double a2_2 = - a2_1;
		double h2_1 = 1.0;
		double h2_2 = h2_1;
		
		double a3_1 = - 0.774596669241483;
		double a3_2 = 0.0;
		double a3_3 = - a3_1;
		double h3_1 = 0.555555555555556;
		double h3_2 = 0.888888888888889;
		double h3_3 = h3_1;
		
		double a4_1 = - 0.861136311594953;
		double a4_2 = - 0.339981043584856;
		double a4_3 = - a4_2;
		double a4_4 = - a4_1;
		double h4_1 = 0.347854845137454;
		double h4_2 = 0.652145154862546;
		double h4_3 = h4_2;
		double h4_4 = h4_1;
		
		double rlt = 0.0;
		if(order == 1) {
			params[paramsStart] = a1_1;
			rlt += h1_1*integrand.apply(params);
		} else if(order == 2) {
			params[paramsStart] = a2_1;
			rlt += h2_1*integrand.apply(params);
			params[paramsStart] = a2_2;
			rlt += h2_2*integrand.apply(params);
		} else if(order == 3) {
			params[paramsStart] = a3_1;
			rlt += h3_1*integrand.apply(params);
			params[paramsStart] = a3_2;
			rlt += h3_2*integrand.apply(params);
			params[paramsStart] = a3_3;
			rlt += h3_3*integrand.apply(params);
		} else if(order == 4) {
			params[paramsStart] = a4_1;
			rlt += h4_1*integrand.apply(params);
			params[paramsStart] = a4_2;
			rlt += h4_2*integrand.apply(params);
			params[paramsStart] = a4_3;
			rlt += h4_3*integrand.apply(params);
			params[paramsStart] = a4_4;
			rlt += h4_4*integrand.apply(params);
		} else if(order == 5) {
			for(int i=0;i<order;i++) {
				params[paramsStart] = a5[i];
				rlt += h5[i]*integrand.apply(params);
			}
		} else {
			System.out.println("ERROR: intOnLinearRefElement() Not supported order = "+order);
		}
		return rlt;
	}

}
