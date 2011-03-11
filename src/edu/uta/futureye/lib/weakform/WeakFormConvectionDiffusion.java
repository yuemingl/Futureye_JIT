package edu.uta.futureye.lib.weakform;

import edu.uta.futureye.core.Element;
import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.basic.SpaceVectorFunction;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.function.intf.VectorFunction;
import edu.uta.futureye.function.operator.FOBasic;
import edu.uta.futureye.function.operator.FOVector;
import edu.uta.futureye.util.Utils;

/**
 * Convection-diffusion equation
 *   \frac{\partial{c}}{\partial{t}} = \Nabla{k*\Nabla{c}} - \mathbf{v}\dot\Nabla{c} + f
 * 
 * Time discrete form:
 *   Let:
 *   \frac{\partial{c}}{\partial{t}} = (c_n+1 - c_n)/Dt
 *   We have,
 *   -Dt*\Nabla{k*\Nabla{c_n+1}} + Dt*\mathbf{v}\dot\Nabla{c_n+1} + c_n+1 = Dt*f + c_n
 * 
 * Weak form:
 *   Let c_n+1 := u
 *   Dt*(k*\Nabla{u},\Nabla{w}) + Dt*( (v1*u_x,w)+(v2*u_y,w)+(v3*u_z,w) ) + b*(u,w) = (Dt*f + c_n,w)
 *   
 * where
 *   c=c(x,y,z,t): particles or energy(e.g. salt density, Heat...) are transferred inside 
 *                 a physical system due to two processes: diffusion and convection
 *   k=k(x,y,z): the diffusion coefficient
 *   \mathbf{v}=(v1,v2,v3)':  the convection velocity vector
 *   f=f(x,y,z): the source term
 *   Dt: the time step size
 *   b: b=1 or b=b(x,y,z), if the equation has the term a*c, where a(x,y,z)=b(x,y,z)-1
 * 
 * Boundary condition
 *   c = c0,                  on \Gamma1 (Dirichlet)
 *   d*c + k*c_n = q,         on \Gamma2 (Robin)
 *
 * The following weak form just gives one step computation of c from c_n to c_n+1. 
 *   
 * @author liuyueming
 *
 */
public class WeakFormConvectionDiffusion extends AbstractScalarWeakForm {
	protected Function g_f = null;
	protected Function g_k = null;
	protected Function g_b = null;
	protected Function g_q = null;
	protected Function g_d = null;
	protected VectorFunction g_v = null;
	protected Function g_cn = null;
	protected double Dt;

	/**
	 * 
	 * @param k
	 * @param b
	 * @param cn: c_n
	 * @param Dt
	 */
	public void setParam(Function k,Function b, Function cn, double Dt) {
		this.g_k = k;
		this.g_b = b;
		this.g_cn = cn;
		this.Dt = Dt;
	}
	
	public void setConvectionVelocity(VectorFunction v) {
		this.g_v = v;
	}
	
	//right hand side function (source term)
	public void setF(Function f) {
		this.g_f = f;
	}
	
	//Robin:  d*u + k*u_n= q (×ÔÈ»±ß½ç£ºd==k, q=0)
	public void setRobin(Function q,Function d) {
		this.g_q = q;
		this.g_d = d;
	}

	@Override
	public Function leftHandSide(Element e, ItemType itemType) {
		 if(itemType==ItemType.Domain)  {
			 //Interplate functions on element e
			Function fk = Utils.interplateFunctionOnElement(g_k,e);
			Function fb = Utils.interplateFunctionOnElement(g_b,e);
			VectorFunction fv = new SpaceVectorFunction(g_v.getDim());
			for(int dim=1;dim<=g_v.getDim();dim++)
				fv.set(dim, Utils.interplateFunctionOnElement(g_v.get(dim),e));
			
			//Dt*(k*\Nabla{u},\Nabla{w}) + Dt*( (v1*u_x,w)+(v2*u_y,w)+(v3*u_z,w) ) + b*(u,w)
			Function integrand = null;
			integrand = fk.X(
							FOVector.Grad(u,u.innerVarNames()).
								dot(
							FOVector.Grad(v,v.innerVarNames()))
						).P(
							fv.dot(FOVector.Grad(u,u.innerVarNames()))
						).X(FC.c(Dt)).P(
							fb.X(u).X(v)
						);
			return integrand;
		}
		else if(itemType==ItemType.Border) {
			if(g_d != null) {
				Element be = e;
				Function fd = Utils.interplateFunctionOnElement(g_d, be);
				Function borderIntegrand = FOBasic.Mult(FOBasic.Mult(fd, u), v);
				return borderIntegrand;
			}
		}
		return null;
	}

	@Override
	public Function rightHandSide(Element e, ItemType itemType) {
		if(itemType==ItemType.Domain)  {
			//(Dt*f + c_n,w)
			Function ff = Utils.interplateFunctionOnElement(g_f, e);
			Function fcn = Utils.interplateFunctionOnElement(g_cn, e);
			Function integrand = FOBasic.Mult(ff.X(FC.c(Dt)).P(fcn),v);
			return integrand;
		} else if(itemType==ItemType.Border) {
			Element be = e;
			Function fq = Utils.interplateFunctionOnElement(g_q, be);
			Function borderIntegrand = FOBasic.Mult(fq, v);
			return borderIntegrand;
		}
		return null;		
	}
}
