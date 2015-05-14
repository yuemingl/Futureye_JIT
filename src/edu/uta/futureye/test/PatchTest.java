package edu.uta.futureye.test;

import edu.uta.futureye.bytecode.CompiledFunc;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.function.basic.FX;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.function.operator.FOIntegrate;
import edu.uta.futureye.lib.shapefun.SFLinearLocal2DRS;
import edu.uta.futureye.util.container.NodeList;
import static edu.uta.futureye.function.FMath.*;

public class PatchTest {

	public static void main(String[] args) {
		NodeList nodes = new NodeList();
		/**
		 * |\
		 * | \
		 * ----
		 */
		nodes.add(new Node(1, 0.0,0.0));
		nodes.add(new Node(2, 0.2,0.0));
		nodes.add(new Node(3, 0.0,0.2));
		Element e = new Element(nodes);
		String[] argsOrder = new String[]{"x1","x2","x3","y1","y2","y3","r","s","t"};
		FX x1 = new FX("x1");
		FX x2 = new FX("x2");
		FX x3 = new FX("x3");
		FX y1 = new FX("y1");
		FX y2 = new FX("y2");
		FX y3 = new FX("y3");
		MathFunc fx = x1*r + x2*s + x3*t;
		MathFunc fy = y1*r + y2*s + y3*t;
		MathFunc f = fx + fy;
		CompiledFunc cf = f.compile(argsOrder);
		double[] params = new double[9];
		double[] coords = e.getNodeCoords();
		System.arraycopy(coords, 0, params, 0, coords.length);
		
		MathFunc fx2 = coords[0]*r + coords[1]*s + coords[2]*t;
		MathFunc fy2 = coords[3]*r + coords[4]*s + coords[5]*t;
		MathFunc f2 = fx2 + fy2;
		
		System.out.println(intOnTriangleRefElement(cf, params, coords.length, 2));
		System.out.println(FOIntegrate.intOnTriangleRefElement(f2, 2));
		
//		SFLinearLocal2D[] shapeFun = new SFLinearLocal2D[3];
//		shapeFun[0] = new SFLinearLocal2D(1);
//		shapeFun[1] = new SFLinearLocal2D(2);
//		shapeFun[2] = new SFLinearLocal2D(3);
		SFLinearLocal2DRS[] shapeFun = new SFLinearLocal2DRS[3];
		shapeFun[0] = new SFLinearLocal2DRS(1);
		shapeFun[1] = new SFLinearLocal2DRS(2);
		shapeFun[2] = new SFLinearLocal2DRS(3);
	}
	
	public static double intOnTriangleRefElement(CompiledFunc integrand, double[] params, int endPos, int order) {
		double rlt = 0.0;
		if(order == 2) {
			params[endPos] = 0.333333333333333;
			params[endPos+1] = 0.333333333333333;
			params[endPos+2] = 0.333333333333333;
			rlt = integrand.apply(params);
		} else if(order == 3) {
			params[endPos] = 0.5; params[endPos+1] = 0.5; params[endPos+2] = 0.0; 
			double pv1 = integrand.apply(params);
			params[endPos] = 0.0; params[endPos+1] = 0.5; params[endPos+2] = 0.5; 
			double pv2 = integrand.apply(params);
			params[endPos] = 0.5; params[endPos+1] = 0.0; params[endPos+2] = 0.5; 
			double pv3 = integrand.apply(params);
			rlt = 0.333333333333333 * (pv1+pv2+pv3);
		}
		return rlt;
//		} else if(order == 4) {
//			double w123 = 25.0/48.0;
//			double w4 = -27.0/48.0;
//			
//			Variable v1 = new Variable();
//			Variable v2 = new Variable();
//			Variable v3 = new Variable();
//			Variable v4 = new Variable();
//			
//			v1.set(VN.r, 0.6); v1.set(VN.s, 0.2); v1.set(VN.t, 0.2);
//			v2.set(VN.r, 0.2); v2.set(VN.s, 0.6); v2.set(VN.t, 0.2);
//			v3.set(VN.r, 0.2); v3.set(VN.s, 0.2); v3.set(VN.t, 0.6);
//			v4.set(VN.r, 0.333333333333333); v4.set(VN.s, 0.333333333333333); v4.set(VN.t, 0.333333333333333);
//			
//			rlt = w123 * (
//					integrand.apply(v1)+
//					integrand.apply(v2)+
//					integrand.apply(v3)
//					) + w4*integrand.apply(v4);
//
//		} else if(order == 5) {
//			Variable v = new Variable();
//			for(int i=0;i<7;i++) {
//				v.set(VN.r, triR[i]); 
//				v.set(VN.s, triS[i]); 
//				v.set(VN.t, 1.0-triR[i]-triS[i]);
//				rlt += triW[i]*integrand.apply(v);
//			}
//		}
//		return 0.5*rlt; //??? 0.5 ???
	}
}
