package edu.uta.futureye.application;

import edu.uta.futureye.core.Element;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.function.operator.FMath;
import edu.uta.futureye.lib.weakform.AbstractScalarWeakForm;
import edu.uta.futureye.util.Utils;

/**
 * Solve Weak form
 *  (k*\nabla{u}, \nabla{v}) + (u*\vec{b}, \nabla{v}) + (c*u, v) = (f, v)
 * =>
 *  (k*u_x, v_x) + (k*u_y, v_y) + (b1*v_x, u) + (b2*v_y, u) + (c*u, v)= (f, v)
 * 
 * @author liuyueming
 *
 */
public class WeakFormGCMDual extends AbstractScalarWeakForm {
	protected MathFunc g_f = null;
	
	protected MathFunc g_k = null;
	protected MathFunc g_c = null;
	protected MathFunc g_b1 = null;
	protected MathFunc g_b2 = null;
	
	protected MathFunc g_q = null;
	protected MathFunc g_d = null;

	public void setF(MathFunc f) {
		this.g_f = f;
	}
	
	public void setParam(MathFunc k,MathFunc c,MathFunc b1,MathFunc b2) {
		this.g_k = k;
		this.g_c = c;
		this.g_b1 = b1;
		this.g_b2 = b2;
	}
	
	//Robin:  d*u + k*u_n= q (自然边界：d==k, q=0)
	public void setRobin(MathFunc q,MathFunc d) {
		this.g_q = q;
		this.g_d = d;
	}	

	@Override
	public MathFunc leftHandSide(Element e, ItemType itemType) {
		if(itemType==ItemType.Domain)  {
			//Integrand part of Weak Form on element e
			MathFunc integrand = null;

			MathFunc fk = Utils.interpolateOnElement(g_k,e);
			MathFunc fc = Utils.interpolateOnElement(g_c,e);
			MathFunc fb1 = Utils.interpolateOnElement(g_b1,e);
			MathFunc fb2 = Utils.interpolateOnElement(g_b2,e);
			
			integrand = FMath.sum(
						fk.M(
								u.diff("x").M(v.diff("x")).A(
								u.diff("y").M(v.diff("y"))
						)),
						fb1.M(v.diff("x").M(u)),
						fb2.M(v.diff("y").M(u)),
						fc.M(u.M(v))
					);
			return integrand;
		}
		else if(itemType==ItemType.Border) {
			if(g_d != null) {
				Element be = e;
				MathFunc fd = Utils.interpolateOnElement(g_d, be);
				MathFunc borderIntegrand = fd.M(u.M(v));
				return borderIntegrand;
			}
		}
		return null;
	}

	@Override
	public MathFunc rightHandSide(Element e, ItemType itemType) {
		if(itemType==ItemType.Domain)  {
			MathFunc ff = Utils.interpolateOnElement(g_f, e);
			MathFunc integrand = ff.M(v);
			return integrand;
		} else if(itemType==ItemType.Border) {
			Element be = e;
			MathFunc fq = Utils.interpolateOnElement(g_q, be);
			MathFunc borderIntegrand = fq.M(v);
			return borderIntegrand;
		} 
		return null;
	}
}
