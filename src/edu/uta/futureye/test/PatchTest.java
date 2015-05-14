package edu.uta.futureye.test;

import edu.uta.futureye.bytecode.CompiledFunc;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.function.basic.FX;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.function.operator.FOIntegrate;
import edu.uta.futureye.lib.shapefun.SFLinearLocal2D;
import edu.uta.futureye.lib.shapefun.SFLinearLocal2DRS;
import edu.uta.futureye.util.container.NodeList;
import static edu.uta.futureye.function.FMath.*;

public class PatchTest {

	public static void shapeFuncTest1() {
		SFLinearLocal2D[] shapeFun = new SFLinearLocal2D[3];
		shapeFun[0] = new SFLinearLocal2D(1);
		shapeFun[1] = new SFLinearLocal2D(2);
		shapeFun[2] = new SFLinearLocal2D(3);
		System.out.println(shapeFun[0]);
		System.out.println(shapeFun[1]);
		System.out.println(shapeFun[2]);
	}
	public static void shapeFuncTest2() {
		SFLinearLocal2DRS[] shapeFun = new SFLinearLocal2DRS[3];
		shapeFun[0] = new SFLinearLocal2DRS(1);
		shapeFun[1] = new SFLinearLocal2DRS(2);
		shapeFun[2] = new SFLinearLocal2DRS(3);
		System.out.println(shapeFun[0]);
		System.out.println(shapeFun[1]);
		System.out.println(shapeFun[2]);
	}
	
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
		
		for(int order=2; order<=5; order++) {
		System.out.println(intOnTriangleRefElement(cf, params, coords.length, order));
		System.out.println(FOIntegrate.intOnTriangleRefElement(f2, order));
		}
	}
	
	public static double intOnTriangleRefElement(CompiledFunc integrand, double[] params, int endPos, int order) {
		double rlt = 0.0;
		if(order == 2) {
			params[endPos] = 0.333333333333333;
			params[endPos+1] = 0.333333333333333;
			params[endPos+2] = 0.333333333333333;
			rlt = 0.5*integrand.apply(params);
		} else if(order == 3) {
			params[endPos] = 0.5; params[endPos+1] = 0.5; params[endPos+2] = 0.0; 
			double pv1 = integrand.apply(params);
			params[endPos] = 0.0; params[endPos+1] = 0.5; params[endPos+2] = 0.5; 
			double pv2 = integrand.apply(params);
			params[endPos] = 0.5; params[endPos+1] = 0.0; params[endPos+2] = 0.5; 
			double pv3 = integrand.apply(params);
			rlt = 0.5*0.333333333333333*(pv1 + pv2 + pv3);
		} else if(order == 4) {
			double w123 = 25.0/48.0;
			double w4 = -27.0/48.0;
			
			params[endPos] = 0.6; params[endPos+1] = 0.2; params[endPos+2] = 0.2; 
			double pv1 = integrand.apply(params);
			params[endPos] = 0.2; params[endPos+1] = 0.6; params[endPos+2] = 0.2; 
			double pv2 = integrand.apply(params);
			params[endPos] = 0.2; params[endPos+1] = 0.2; params[endPos+2] = 0.6; 
			double pv3 = integrand.apply(params);
			params[endPos] = 0.333333333333333; params[endPos+1] = 0.333333333333333; params[endPos+2] = 0.333333333333333; 
			double pv4 = 0.5*integrand.apply(params);
			
			rlt = 0.5*w123*(pv1 + pv2 + pv3) + w4*pv4;
		} else if(order == 5) {
			for(int i=0;i<7;i++) {
				params[endPos]   = triR[i]; 
				params[endPos+1] = triS[i]; 
				params[endPos+2] = 1.0-triR[i]-triS[i];
				rlt += triW[i]*integrand.apply(params);
			}
		}
		return rlt;
	}
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
}
