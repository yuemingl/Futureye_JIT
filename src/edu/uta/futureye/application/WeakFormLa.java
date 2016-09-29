package edu.uta.futureye.application;

import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.core.geometry.Point;
import edu.uta.futureye.function.FMath;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.FXY;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.lib.weakform.AbstractScalarWeakForm;
import edu.uta.futureye.util.Utils;

/**
 * Solve u, weak form:
 * 
 *  (u, v) = (f + k1*\sum{\nabla{b_i}\cdot\nabla{c_i}} + k2*\sum{b_i*c_i}, v), for all v
 *  =>
 *  (u, v) = (f + k1*\sum{(b_ix*c_ix+b_iy*c_iy)} + k2*\sum{b_i*c_i}, v), for all v
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
public class WeakFormLa extends AbstractScalarWeakForm {
	protected MathFunc g_f = null;
	
	protected MathFunc g_k1 = null;
	protected MathFunc g_k2 = null;
	protected MathFunc[] g_b = null;
	protected MathFunc[] g_c = null;
	
	protected MathFunc g_q = null;
	protected MathFunc g_d = null;

	public void setF(MathFunc f,
			MathFunc k1,MathFunc k2,
			MathFunc[] b,MathFunc[] c) {
		this.g_f = f;
		this.g_c = c;
		this.g_b = b;
		this.g_k1 = k1;		
		this.g_k2 = k2;		
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
			integrand = u.M(v);
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
			MathFunc fk1 = Utils.interpolateOnElement(g_k1,e);
			MathFunc fk2 = Utils.interpolateOnElement(g_k2,e);
			int NF = g_b.length;
			MathFunc[] fb = new MathFunc[NF];
			MathFunc[] fc = new MathFunc[NF];
			MathFunc[] fbx = new MathFunc[NF];
			MathFunc[] fby = new MathFunc[NF];
			MathFunc[] fcx = new MathFunc[NF];
			MathFunc[] fcy = new MathFunc[NF];
			MathFunc[] fbxMfcx = new MathFunc[NF];
			MathFunc[] fbyMfcy = new MathFunc[NF];
			MathFunc[] fbMfc = new MathFunc[NF];

			for(int k=0;k<NF;k++) {
				fb[k] = Utils.interpolateOnElement(g_b[k],e);
				fc[k] = Utils.interpolateOnElement(g_c[k],e);

				int N = e.nodes.size();
				double[] fbv = new double[N];
				double[] fcv = new double[N];
				for(int i=1;i<=N;i++) {
					Node node = e.nodes.at(i);
					Variable var = Variable.createFrom(g_b[k], node, node.globalIndex);
					fbv[i-1] = g_b[k].apply(var);
					fcv[i-1] = g_c[k].apply(var);
				}
				//d(a1 + a2*x + a3*y + a4*x*y)/dx
				//d(a1 + a2*x + a3*y + a4*x*y)/dy
				double[] ab = Utils.computeBilinearFunctionCoef(e.nodes.toArray(new Point[0]), fbv);
				double[] ac = Utils.computeBilinearFunctionCoef(e.nodes.toArray(new Point[0]), fcv);
				MathFunc bdx = new FXY(0.0,ab[3],ab[1]);
				MathFunc bdy = new FXY(ab[3],0.0,ab[2]);
				MathFunc cdx = new FXY(0.0,ac[3],ac[1]);
				MathFunc cdy = new FXY(ac[3],0.0,ac[2]);
			
				fbx[k] = Utils.interpolateOnElement(bdx, e);
				fby[k] = Utils.interpolateOnElement(bdy, e);
				fcx[k] = Utils.interpolateOnElement(cdx, e);
				fcy[k] = Utils.interpolateOnElement(cdy, e);
				
				fbxMfcx[k] = fbx[k].M(fcx[k]);
				fbyMfcy[k] = fby[k].M(fcy[k]);
				fbMfc[k] = fb[k].M(fc[k]);
			}
			
			MathFunc integrand = FMath.sum(
						ff,
						fk1.M(FMath.sum(fbxMfcx)),
						fk1.M(FMath.sum(fbyMfcy)),
						fk2.M(FMath.sum(fbMfc))
					).M(v);
			
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
