package edu.uta.futureye.application;

import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.core.geometry.Point;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.FXY;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.function.operator.FMath;
import edu.uta.futureye.lib.weakform.AbstractScalarWeakForm;
import edu.uta.futureye.util.Utils;

/**
 * Solve u, weak form:
 * 
 *  (k*u*\nabla{b}, \nabla{v}) + (c*u, v)= (f, v), for all v
 *  =>
 *  (k*u*b_x, v_x) + (k*u*b_y, v_y) + (c*u, v)= (f, v), for all v
 *  
 * where
 *   k=k(x,y)
 *   b=b(x,y)
 *   \nabla{b} = (b_x, b_y)
 *   c=c(x,y)
 *   f=f(x,y)
 *    
 * @author liuyueming
 *
 */
public class WeakFormC extends AbstractScalarWeakForm {
	protected Function g_f = null;
	
	protected Function g_k = null;
	protected Function g_c = null;
	protected Function g_b = null;
	
	protected Function g_q = null;
	protected Function g_d = null;

	public void setF(Function f) {
		this.g_f = f;
	}
	
	public void setParam(Function k,Function c,Function b) {
		this.g_k = k;
		this.g_c = c;
		this.g_b = b;
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

			Function fk = Utils.interpolateFunctionOnElement(g_k,e);
			Function fc = Utils.interpolateFunctionOnElement(g_c,e);
			
			int N = e.nodes.size();
			double[] f = new double[N];
			for(int i=1;i<=N;i++) {
				Node node = e.nodes.at(i);
				Variable var = Variable.createFrom(g_b, node, node.globalIndex);
				f[i-1] = g_b.value(var);
			}
			double[] a = Utils.computeBilinearFunctionCoef(e.nodes.toArray(new Point[0]), f);
			//d(a1 + a2*x + a3*y + a4*x*y)/dx
			Function dx = new FXY(0.0,a[3],a[1]);
			//d(a1 + a2*x + a3*y + a4*x*y)/dy
			Function dy = new FXY(a[3],0.0,a[2]);
			
			Function fbx = Utils.interpolateFunctionOnElement(dx, e);
			Function fby = Utils.interpolateFunctionOnElement(dy, e);
			
			integrand = FMath.sum(
						fk.M(fbx.M(u.M(v._d("x")))),
						fk.M(fby.M(u.M(v._d("y")))),
						fc.M(u.M(v))
					);
			return integrand;
		}
		else if(itemType==ItemType.Border) {
			if(g_d != null) {
				Element be = e;
				Function fd = Utils.interpolateFunctionOnElement(g_d, be);
				Function borderIntegrand = fd.M(u.M(v));
				return borderIntegrand;
			}
		}
		return null;
	}

	@Override
	public Function rightHandSide(Element e, ItemType itemType) {
		if(itemType==ItemType.Domain)  {
			Function ff = Utils.interpolateFunctionOnElement(g_f, e);
			Function integrand = ff.M(v);
			return integrand;
		} else if(itemType==ItemType.Border) {
			Element be = e;
			Function fq = Utils.interpolateFunctionOnElement(g_q, be);
			Function borderIntegrand = fq.M(v);
			return borderIntegrand;
		} 
		return null;
	}
}
