package edu.uta.futureye.lib.weakform;

import edu.uta.futureye.core.Element;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.function.operator.FMath;
import edu.uta.futureye.util.Utils;

/**
 * Solve
 *   -k*Laplace(u) + c*u = f, in \Omega
 *   u = u0,                  on \Gamma1
 *   d*u + k*u_n = q,         on \Gamma2
 *=>
 *   A(u, v) = (f, v)
 * 
 * where
 *   A(u, v) = (k*Grad{u}, Grad{v}) - (q-d*u,v)_\Gamma2 + (c*u, v)
 *
 *   \Gamma1: Dirichlet boundary of \Omega
 *   \Gamma2: Neumann(Robin) boundary of \Omega
 *   u_n: \frac{\pratial{u}}{\partial{n}}
 *   n: unit norm vector of \Omega
 *   k = k(x,y)
 *   c = c(x,y)
 *   d = d(x,y)
 *   q = q(x,y)
 *   
 * @author liuyueming
 *
 */
public class WeakFormLaplace extends AbstractScalarWeakForm {
	protected Function g_f = null;
	protected Function g_k = null;
	protected Function g_c = null;
	protected Function g_q = null;
	protected Function g_d = null;

	//right hand side function (source term)
	public void setF(Function f) {
		this.g_f = f;
	}
	
	//Robin: k*u_n + d*u = q (自然边界：d==c)
	public void setParam(Function k,Function c,Function q,Function d) {
		this.g_k = k;
		this.g_c = c;
		this.g_q = q;
		this.g_d = d;
	}

	@Override
	public Function leftHandSide(Element e, ItemType itemType) {
		if(itemType==ItemType.Domain)  {
			//Integrand part of Weak Form on element e
			Function integrand = null;
			if(g_k == null) {
				integrand =  
					FMath.grad(u,u.innerVarNames()).
								  dot( 
					FMath.grad(v,v.innerVarNames()) 
								  );
			} else {
				
				Function fk = Utils.interplateFunctionOnElement(g_k,e);
				Function fc = Utils.interplateFunctionOnElement(g_c,e);
				integrand = fk.M(
									FMath.grad(u,u.innerVarNames()).
									dot(
									FMath.grad(v,v.innerVarNames())))
							.A(
									fc.M(u.M(v))
							);
			}
			return integrand;
		}
		else if(itemType==ItemType.Border) {
			if(g_d != null) {
				Element be = e;
				Function fd = Utils.interplateFunctionOnElement(g_d, be);
				Function borderIntegrand = fd.M(u.M(v));
				return borderIntegrand;
			}
		}
		return null;
	}

	@Override
	public Function rightHandSide(Element e, ItemType itemType) {
		if(itemType==ItemType.Domain)  {
			Function ff = Utils.interplateFunctionOnElement(g_f, e);
			Function integrand = ff.M(v);
			return integrand;
		} else if(itemType==ItemType.Border) {
			Element be = e;
			Function fq = Utils.interplateFunctionOnElement(g_q, be);
			Function borderIntegrand = fq.M(v);
			return borderIntegrand;
		}
		return null;		
	}
}
