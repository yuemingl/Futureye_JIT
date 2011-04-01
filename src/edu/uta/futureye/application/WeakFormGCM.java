package edu.uta.futureye.application;

import edu.uta.futureye.core.Element;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.function.operator.FMath;
import edu.uta.futureye.lib.weakform.AbstractScalarWeakForm;
import edu.uta.futureye.util.Utils;

/**
 * Solve
 *  -\Delta{u} + \mathbf{b}\dot\Nabla{u} + c*u = f
 * 
 * Weak form
 *  (k*u_x, v_x) + (k*u_y, v_y) + (b1*u_x, v) + (b2*u_y, v) + (c*u, v)= (f, v)
 * 
 * @author liuyueming
 *
 */
public class WeakFormGCM extends AbstractScalarWeakForm {
	protected Function g_f = null;
	
	protected Function g_k = null;
	protected Function g_c = null;
	protected Function g_b1 = null;
	protected Function g_b2 = null;
	
	protected Function g_q = null;
	protected Function g_d = null;

	public void setF(Function f) {
		this.g_f = f;
	}
	
	public void setParam(Function k,Function c,Function b1,Function b2) {
		this.g_k = k;
		this.g_c = c;
		this.g_b1 = b1;
		this.g_b2 = b2;
	}
	
	//Robin:  d*u + k*u_n= q (自然边界：d==k, q=0)
	public void setRobin(Function q,Function d) {
		this.g_q = q;
		this.g_d = d;
	}	

	@Override
	public Function leftHandSide(Element e, ItemType itemType) {
		if(itemType==ItemType.Domain)  {
			//Integrand part of Weak Form on element e
			Function integrand = null;

			Function fk = Utils.interplateFunctionOnElement(g_k,e);
			Function fc = Utils.interplateFunctionOnElement(g_c,e);
			Function fb1 = Utils.interplateFunctionOnElement(g_b1,e);
			Function fb2 = Utils.interplateFunctionOnElement(g_b2,e);
			
			integrand = FMath.PlusAll(
						FMath.Mult(fk, FMath.Plus(
								FMath.Mult(u._d("x"), v._d("x")),
								FMath.Mult(u._d("y"), v._d("y"))
						)),
						FMath.Mult(fb1,FMath.Mult(u._d("x"), v)),
						FMath.Mult(fb2,FMath.Mult(u._d("y"), v)),
						FMath.Mult(fc, FMath.Mult(u, v))
					);
			return integrand;
		}
		else if(itemType==ItemType.Border) {
			if(g_d != null) {
				Element be = e;
				Function fd = Utils.interplateFunctionOnElement(g_d, be);
				Function borderIntegrand = FMath.Mult(FMath.Mult(fd, u), v);
				return borderIntegrand;
			}
		}
		return null;
	}

	@Override
	public Function rightHandSide(Element e, ItemType itemType) {
		if(itemType==ItemType.Domain)  {
			Function ff = Utils.interplateFunctionOnElement(g_f, e);
			Function integrand = FMath.Mult(ff,v);
			return integrand;
		} else if(itemType==ItemType.Border) {
			Element be = e;
			Function fq = Utils.interplateFunctionOnElement(g_q, be);
			Function borderIntegrand = FMath.Mult(fq, v);
			return borderIntegrand;
		} 
		return null;
	}
}
