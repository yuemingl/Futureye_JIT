package edu.uta.futureye.application;

import edu.uta.futureye.core.Element;
import edu.uta.futureye.function.intf.MathFun;
import edu.uta.futureye.function.operator.FMath;
import edu.uta.futureye.lib.weakform.AbstractScalarWeakForm;
import edu.uta.futureye.util.Utils;

/**
 * Solve
 *  -k\Delta{u} + \mathbf{b}\cdot\nabla{u} + c*u = f
 * 
 * Weak form
 *  (k*u_x, v_x) + (k*u_y, v_y) + (b1*u_x, v) + (b2*u_y, v) + (c*u, v)= (f, v)
 * 
 * @author liuyueming
 *
 */
public class WeakFormGCM extends AbstractScalarWeakForm {
	protected MathFun g_f = null;
	
	protected MathFun g_k = null;
	protected MathFun g_c = null;
	protected MathFun g_b1 = null;
	protected MathFun g_b2 = null;
	
	protected MathFun g_q = null;
	protected MathFun g_d = null;

	public void setF(MathFun f) {
		this.g_f = f;
	}
	
	public void setParam(MathFun k,MathFun c,MathFun b1,MathFun b2) {
		this.g_k = k;
		this.g_c = c;
		this.g_b1 = b1;
		this.g_b2 = b2;
	}
	
	//Robin:  d*u + k*u_n= q (自然边界：d==k, q=0)
	public void setRobin(MathFun q,MathFun d) {
		this.g_q = q;
		this.g_d = d;
	}	

	@Override
	public MathFun leftHandSide(Element e, ItemType itemType) {
		if(itemType==ItemType.Domain)  {
			//Integrand part of Weak Form on element e
			MathFun integrand = null;

			MathFun fk = Utils.interpolateOnElement(g_k,e);
			MathFun fc = Utils.interpolateOnElement(g_c,e);
			MathFun fb1 = Utils.interpolateOnElement(g_b1,e);
			MathFun fb2 = Utils.interpolateOnElement(g_b2,e);
			
			integrand = FMath.sum(
						fk.M(
								u._d("x").M(v._d("x")).A(
								u._d("y").M(v._d("y"))
						)),
						fb1.M(u._d("x").M(v)),
						fb2.M(u._d("y").M(v)),
						fc.M(u.M(v))
					);
			return integrand;
		}
		else if(itemType==ItemType.Border) {
			if(g_d != null) {
				Element be = e;
				MathFun fd = Utils.interpolateOnElement(g_d, be);
				MathFun borderIntegrand = fd.M(u.M(v));
				return borderIntegrand;
			}
		}
		return null;
	}

	@Override
	public MathFun rightHandSide(Element e, ItemType itemType) {
		if(itemType==ItemType.Domain)  {
			MathFun ff = Utils.interpolateOnElement(g_f, e);
			MathFun integrand = ff.M(v);
			return integrand;
		} else if(itemType==ItemType.Border) {
			Element be = e;
			MathFun fq = Utils.interpolateOnElement(g_q, be);
			MathFun borderIntegrand = fq.M(v);
			return borderIntegrand;
		} 
		return null;
	}
}
