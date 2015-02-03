package edu.uta.futureye.lib.weakform;

import edu.uta.futureye.core.Element;
import edu.uta.futureye.function.intf.MathFun;
import edu.uta.futureye.util.Utils;
import static edu.uta.futureye.function.operator.FMath.*;

/**
 * <blockquote><pre>
 * Solve the following problem(2D or 3D):
 *   -k*Laplace(u) + c*u = f, in \Omega
 *   u = u0,                  on \Gamma1
 *   d*u + k*u_n = g,         on \Gamma2
 *   
 *=>Weak formulation:
 *   A(u, v) = (f, v)
 *   
 * where
 *   A(u, v) = (k*Grad{u}, Grad{v}) - (g-d*u,v)_\Gamma2 + (c*u, v)
 *   \Gamma1: Dirichlet boundary of \Omega
 *   \Gamma2: Neumann(Robin) boundary of \Omega
 *   u_n: \frac{\pratial{u}}{\partial{n}}
 *   \vec{n}: unit norm vector of \Omega
 *   k = k(\vec{x})
 *   c = c(\vec{x})
 *   d = d(\vec{x})
 *   g = q(\vec{x})
 *
 * Remark:
 * *For nature boundary condition (自然边界条件)
 *   k*u_n + d*u = 0, set d=k
 * =>
 *   u_n + u = 0 
 * </blockquote></pre>  
 *   
 * @author liuyueming
 *
 */
public class WeakFormLaplace extends AbstractScalarWeakForm {
	protected MathFun g_f = null;
	protected MathFun g_k = null;
	protected MathFun g_c = null;
	protected MathFun g_g = null;
	protected MathFun g_d = null;

	//right hand side function (source term)
	public void setF(MathFun f) {
		this.g_f = f;
	}
	
	//Robin: d*u +  k*u_n = g
	public void setParam(MathFun k,MathFun c,MathFun g,MathFun d) {
		this.g_k = k;
		this.g_c = c;
		this.g_g = g;
		this.g_d = d;
	}

	@Override
	public MathFun leftHandSide(Element e, ItemType itemType) {
		if(itemType==ItemType.Domain)  {
			//Integrand part of Weak Form on element e
			MathFun integrand = null;
			if(g_k == null) {
				integrand = grad(u,u.innerVarNames()).dot(
							grad(v,v.innerVarNames()));
			} else {
				MathFun fk = Utils.interpolateOnElement(g_k,e);
				MathFun fc = Utils.interpolateOnElement(g_c,e);
				integrand = fk.M(
					grad(u,u.innerVarNames()).dot(
					grad(v,v.innerVarNames()))).A(
					fc.M(u.M(v)));
			}
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
			if(g_g != null) {
				Element be = e;
				MathFun fq = Utils.interpolateOnElement(g_g, be);
				MathFun borderIntegrand = fq.M(v);
				return borderIntegrand;
			}
		}
		return null;		
	}
}
