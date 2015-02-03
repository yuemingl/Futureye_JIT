package edu.uta.futureye.application;

import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.core.geometry.Point;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.FXY;
import edu.uta.futureye.function.intf.MathFun;
import edu.uta.futureye.function.operator.FMath;
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
	protected MathFun g_f = null;
	
	protected MathFun g_k1 = null;
	protected MathFun g_k2 = null;
	protected MathFun[] g_b = null;
	protected MathFun[] g_c = null;
	
	protected MathFun g_q = null;
	protected MathFun g_d = null;

	public void setF(MathFun f,
			MathFun k1,MathFun k2,
			MathFun[] b,MathFun[] c) {
		this.g_f = f;
		this.g_c = c;
		this.g_b = b;
		this.g_k1 = k1;		
		this.g_k2 = k2;		
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
			integrand = u.M(v);
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
			MathFun fk1 = Utils.interpolateOnElement(g_k1,e);
			MathFun fk2 = Utils.interpolateOnElement(g_k2,e);
			int NF = g_b.length;
			MathFun[] fb = new MathFun[NF];
			MathFun[] fc = new MathFun[NF];
			MathFun[] fbx = new MathFun[NF];
			MathFun[] fby = new MathFun[NF];
			MathFun[] fcx = new MathFun[NF];
			MathFun[] fcy = new MathFun[NF];
			MathFun[] fbxMfcx = new MathFun[NF];
			MathFun[] fbyMfcy = new MathFun[NF];
			MathFun[] fbMfc = new MathFun[NF];

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
				MathFun bdx = new FXY(0.0,ab[3],ab[1]);
				MathFun bdy = new FXY(ab[3],0.0,ab[2]);
				MathFun cdx = new FXY(0.0,ac[3],ac[1]);
				MathFun cdy = new FXY(ac[3],0.0,ac[2]);
			
				fbx[k] = Utils.interpolateOnElement(bdx, e);
				fby[k] = Utils.interpolateOnElement(bdy, e);
				fcx[k] = Utils.interpolateOnElement(cdx, e);
				fcy[k] = Utils.interpolateOnElement(cdy, e);
				
				fbxMfcx[k] = fbx[k].M(fcx[k]);
				fbyMfcy[k] = fby[k].M(fcy[k]);
				fbMfc[k] = fb[k].M(fc[k]);
			}
			
			MathFun integrand = FMath.sum(
						ff,
						fk1.M(FMath.sum(fbxMfcx)),
						fk1.M(FMath.sum(fbyMfcy)),
						fk2.M(FMath.sum(fbMfc))
					).M(v);
			
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
