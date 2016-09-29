package edu.uta.futureye.lib.weakform;

import edu.uta.futureye.core.Element;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.function.intf.VectorFunction;
import edu.uta.futureye.util.FutureyeException;
import edu.uta.futureye.util.Utils;

public class WeakFormElasticIsoPlaneStress2D extends AbstractVectorWeakForm {
	double E = 1000; //Young's modulus
	double gamma = 0.25; //Poisson ratio
	VectorFunction g_b; //Body force
	VectorFunction g_t; //Distributed external loading on boundary
	
	@Override
	public MathFunc leftHandSide(Element e, ItemType itemType) {
		if(itemType==ItemType.Domain)  {
			MathFunc u1 = u.get(1);
			MathFunc u2 = u.get(2);
			MathFunc v1 = v.get(1);
			MathFunc v2 = v.get(2);
			MathFunc u1x = u1.diff("x");
			MathFunc u1y = u1.diff("y");
			MathFunc u2x = u2.diff("x");
			MathFunc u2y = u2.diff("y");
			MathFunc v1x = v1.diff("x");
			MathFunc v1y = v1.diff("y");
			MathFunc v2x = v2.diff("x");
			MathFunc v2y = v2.diff("y");
			
			double coef1 = E/(1-gamma*gamma);
			double coef2 = (1-gamma)/2.0;
			MathFunc integrand = null;
			
//			integrand = 
//			u1x.M(u1x.A(u2y.M(gamma)))  .A(
//			u2y.M(u1x.M(gamma).A(u2y)) ).A(
//			u1y.A(u2x).M(u1y.A(u2x)).M(coef2)								 
//			).M(coef1);
			
			//自由度双循环的表达规则
			//e.g. \vec{u}=(u v)
			//u_x*u_x = u1x*v1x
			//u_y*u_y = u1y*v1y
			//v_x*v_x = u2x*v2x
			//v_y*v_y = u2y*v2y
			//u*v = u1*v2 or u2*v1
			//u_x*v_y = u1x*v2y or u2x*v1y
			//u_y*v_x = u1y*v2x or u2y*v1x
			integrand = 
			u1x.M(v1x).A(
			u2y.M(v2y)).A(
			u1x.M(v2y).M(gamma)
			).A(
				u1y.M(v1y).A(
				u2x.M(v2x).A(
				v1y.M(u2x).M(2)
				)).M(coef2)
			).M(coef1); 

			
			return integrand;
		}
		return null;
	}

	@Override
	public MathFunc rightHandSide(Element e, ItemType itemType) {
		if(itemType==ItemType.Domain)  {
			VectorFunction fb = Utils.interpolateOnElement(g_b, e);
			return v.dot(fb);
		} else if(itemType==ItemType.Border) {
			Element be = e;
			VectorFunction ft = Utils.interpolateOnElement(g_t, be);
			MathFunc borderIntegrand = v.dot(ft);
			return borderIntegrand;
		}
		return null;
	}

	/**
	 * Body force
	 * @param f
	 */
	public void setF(VectorFunction b, VectorFunction t) {
		this.g_b = b;
		this.g_t = t;
	}
	
	/**
	 * 
	 * @param E
	 * @param gamma
	 */
	public void setParam(double E, double gamma) {
		this.E = E;
		this.gamma = gamma;
	}
	
	public boolean isVVFComponentCoupled(int nComponent1, int nComponent2) {
		return true;
	}
	
}
