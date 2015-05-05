package edu.uta.futureye.application;

import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.core.geometry.Point;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.FXY;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.function.operator.FMath;
import edu.uta.futureye.lib.weakform.AbstractScalarWeakForm;
import edu.uta.futureye.util.Utils;

/**
 * Solve u, weak form:
 * 
 *  (k*\nabla{b}*\nabla{u}, v) + (c*u, v)= (f, v), for all v
 *  =>
 *  (k*b_x*u_x, v) + (k*b_y*u_y, v) + (c*u, v)= (f, v), for all v
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
public class WeakFormCT extends AbstractScalarWeakForm {
	protected MathFunc g_f = null;
	
	protected MathFunc g_k = null;
	protected MathFunc g_c = null;
	protected MathFunc g_b = null;
	
	protected MathFunc g_q = null;
	protected MathFunc g_d = null;

	public void setF(MathFunc f) {
		this.g_f = f;
	}
	
	public void setParam(MathFunc k,MathFunc c,MathFunc b) {
		this.g_k = k;
		this.g_c = c;
		this.g_b = b;
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
			
			int N = e.nodes.size();
			double[] f = new double[N];
			for(int i=1;i<=N;i++) {
				Node node = e.nodes.at(i);
				Variable var = Variable.createFrom(g_b, node, node.globalIndex);
				f[i-1] = g_b.apply(var);
			}
			double[] a = Utils.computeBilinearFunctionCoef(e.nodes.toArray(new Point[0]), f);
			//d(a1 + a2*x + a3*y + a4*x*y)/dx
			MathFunc dx = new FXY(0.0,a[3],a[1]);
			//d(a1 + a2*x + a3*y + a4*x*y)/dy
			MathFunc dy = new FXY(a[3],0.0,a[2]);
			
			MathFunc fbx = Utils.interpolateOnElement(dx, e);
			MathFunc fby = Utils.interpolateOnElement(dy, e);
			
			integrand = FMath.sum(
						fk.M(fbx.M(u._d("x").M(v))),
						fk.M(fby.M(u._d("y").M(v))),
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
