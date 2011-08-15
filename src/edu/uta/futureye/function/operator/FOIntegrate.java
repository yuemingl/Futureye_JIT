package edu.uta.futureye.function.operator;

import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.intf.Function;

public class FOIntegrate{
	/**
	 * 在参考三角形单元上积分
	 * @param integrand
	 * @param degree
	 * @return Function
	 */	
	public static double intOnTriangleRefElement(Function integrand, int degree) {
		double dOneThird = 1.0D/3.0D;
		double dOneTwo = 0.5D;

		double rlt = 0.0;
		if(degree == 2) {
			Variable v = new Variable();
			v.set("r", dOneThird);
			v.set("s", dOneThird);
			v.set("t", dOneThird);
			rlt = integrand.value(v);
		} else if(degree == 3) {
			
			Variable v1 = new Variable();
			Variable v2 = new Variable();
			Variable v3 = new Variable();
			
			v1.set("r", dOneTwo); v1.set("s", dOneTwo); v1.set("t", 0.0);
			v2.set("r", 0.0); v2.set("s", dOneTwo); v2.set("t", dOneTwo);
			v3.set("r", dOneTwo); v3.set("s", 0.0); v3.set("t", dOneTwo);
			double pv1 = integrand.value(v1);
			double pv2 = integrand.value(v2);
			double pv3 = integrand.value(v3);
			
			rlt = dOneThird * (pv1+pv2+pv3);
		} else if(degree == 4) {
			double w123 = 25.0/48.0;
			double w4 = -27.0/48.0;
			
			Variable v1 = new Variable();
			Variable v2 = new Variable();
			Variable v3 = new Variable();
			Variable v4 = new Variable();
			
			v1.set("r", 0.6); v1.set("s", 0.2); v1.set("t", 0.2);
			v2.set("r", 0.2); v2.set("s", 0.6); v2.set("t", 0.2);
			v3.set("r", 0.2); v3.set("s", 0.2); v3.set("t", 0.6);
			v4.set("r", dOneThird); v4.set("s", dOneThird); v4.set("t", dOneThird);
			
			rlt = w123 * (
					integrand.value(v1)+
					integrand.value(v2)+
					integrand.value(v3)
					) + w4*integrand.value(v4);
			
		} else if(degree == 5) {
	        double[] w = {
	        	0.06296959,
		        0.06619708,
		        0.06296959,
		        0.06619708,
		        0.06296959,
		        0.06619708,
		        0.11250000
	        };
		    double[]x = {
		        0.10128651,
		        0.47014206,
		        0.79742699,
		        0.47014206,
		        0.10128651,
		        0.05971587,
		        1.0 / 3.0
	        };
		    double[]y = {
				0.10128651,
				0.05971587,
				0.10128651,
				0.47014206,
				0.79742699,
				0.47014206,
				1.0 / 3.0
		    };
		    
			Variable v = new Variable();
			for(int i=0;i<7;i++) {
				v.set("r", x[i]); 
				v.set("s", y[i]); 
				v.set("t", 1.0-x[i]-y[i]);
				rlt += w[i]*integrand.value(v);
			}	
		}
		//???  0.5 ???
		return 0.5*rlt;
	}
	
	
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
	
	/**
	 * 一维参考单元[-1,1]上的数值积分
	 * @param integrand
	 * @param degree
	 * @return
	 */
	public static double intOnLinearRefElement(Function integrand, int degree) {
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
		if(degree == 1) {
			v.set("r", a1_1);
			rlt += h1_1*integrand.value(v);
		} else if(degree == 2) {
			v.set("r", a2_1);
			rlt += h2_1*integrand.value(v);
			v.set("r", a2_2);
			rlt += h2_2*integrand.value(v);
		} else if(degree == 3) {
			v.set("r", a3_1);
			rlt += h3_1*integrand.value(v);
			v.set("r", a3_2);
			rlt += h3_2*integrand.value(v);
			v.set("r", a3_3);
			rlt += h3_3*integrand.value(v);
		} else if(degree == 4) {
			v.set("r", a4_1);
			rlt += h4_1*integrand.value(v);
			v.set("r", a4_2);
			rlt += h4_2*integrand.value(v);
			v.set("r", a4_3);
			rlt += h4_3*integrand.value(v);
			v.set("r", a4_4);
			rlt += h4_4*integrand.value(v);
		} else if(degree == 5) {
			for(int i=0;i<degree;i++) {
				v.set("r", a5[i]);
				rlt += h5[i]*integrand.value(v);
			}
		} else {
			System.out.println("ERROR: intOnLinearRefElement() Not supported degree = "+degree);
		}
		return rlt;
	}
	
	public static double intOnRectangleRefElement(Function integrand, int degree) {
		double a1_1 = 0.0;
		double h1_1 = 4.0;
		double a2 = 0.577350269189626;
		double h2 = 1.0;
		
		Variable v = new Variable();
		double rlt = 0.0;
		if(degree == 1) {
			v.set("r", a1_1);
			v.set("s", a1_1);
			rlt += h1_1*integrand.value(v);
		} else if(degree == 2) {
			v.set("r", a2);
			v.set("s", a2);
			rlt += h2*integrand.value(v);			
			v.set("r", -a2);
			v.set("s", a2);
			rlt += h2*integrand.value(v);			
			v.set("r", a2);
			v.set("s", -a2);
			rlt += h2*integrand.value(v);			
			v.set("r", -a2);
			v.set("s", -a2);
			rlt += h2*integrand.value(v);			
		} else if(degree == 5) {
			for(int i=0;i<degree;i++) {
				for(int j=0;j<degree;j++) {
					v.set("r", a5[i]);
					v.set("s", a5[j]);
					rlt += h5[i]*h5[j]*integrand.value(v);
				}
			}
		} else {
			System.out.println("ERROR: intOnLinearRefElement() Not supported degree = "+degree);
		}
		
		return rlt;
	}		
	
	public static double intOnTetrahedraRefElement(Function integrand, int degree) {
		double a1_1 = 0.25;
		double h1_1 = 1;
		
		double []a2 = {0.585410196624969,0.138196601125011,0.138196601125011,0.138196601125011};
		double h2 = 0.25;
		int [][]M24 = {{0,1,2,3},{1,0,2,3},{1,2,0,3},{1,2,3,0}};
		
		Variable v = new Variable();
		double rlt = 0.0;
		if(degree == 1) {
			v.set("r", a1_1);
			v.set("s", a1_1);
			v.set("t", a1_1);
			v.set("u", a1_1);
			rlt += h1_1*integrand.value(v);
		} else if(degree ==2) {
			for(int i=0;i<M24.length;i++) {
				v.set("r", a2[M24[i][0]]);
				v.set("s", a2[M24[i][1]]);
				v.set("t", a2[M24[i][2]]);
				v.set("u", a2[M24[i][3]]);
				rlt += h2*integrand.value(v);			
			}
		} else {
			System.out.println("ERROR: intOnLinearRefElement() Not supported degree = "+degree);
		}
		
		return rlt;
	}		
}
