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
	protected Function g_f = null;
	
	protected Function g_k1 = null;
	protected Function g_k2 = null;
	protected Function[] g_b = null;
	protected Function[] g_c = null;
	
	protected Function g_q = null;
	protected Function g_d = null;

	public void setF(Function f,
			Function k1,Function k2,
			Function[] b,Function[] c) {
		this.g_f = f;
		this.g_c = c;
		this.g_b = b;
		this.g_k1 = k1;		
		this.g_k2 = k2;		
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
			integrand = u.M(v);
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
			Function fk1 = Utils.interpolateFunctionOnElement(g_k1,e);
			Function fk2 = Utils.interpolateFunctionOnElement(g_k2,e);
			int NF = g_b.length;
			Function[] fb = new Function[NF];
			Function[] fc = new Function[NF];
			Function[] fbx = new Function[NF];
			Function[] fby = new Function[NF];
			Function[] fcx = new Function[NF];
			Function[] fcy = new Function[NF];
			Function[] fbxMfcx = new Function[NF];
			Function[] fbyMfcy = new Function[NF];
			Function[] fbMfc = new Function[NF];

			for(int k=0;k<NF;k++) {
				fb[k] = Utils.interpolateFunctionOnElement(g_b[k],e);
				fc[k] = Utils.interpolateFunctionOnElement(g_c[k],e);

				int N = e.nodes.size();
				double[] fbv = new double[N];
				double[] fcv = new double[N];
				for(int i=1;i<=N;i++) {
					Node node = e.nodes.at(i);
					Variable var = Variable.createFrom(g_b[k], node, node.globalIndex);
					fbv[i-1] = g_b[k].value(var);
					fcv[i-1] = g_c[k].value(var);
				}
				//d(a1 + a2*x + a3*y + a4*x*y)/dx
				//d(a1 + a2*x + a3*y + a4*x*y)/dy
				double[] ab = Utils.computeBilinearFunctionCoef(e.nodes.toArray(new Point[0]), fbv);
				double[] ac = Utils.computeBilinearFunctionCoef(e.nodes.toArray(new Point[0]), fcv);
				Function bdx = new FXY(0.0,ab[3],ab[1]);
				Function bdy = new FXY(ab[3],0.0,ab[2]);
				Function cdx = new FXY(0.0,ac[3],ac[1]);
				Function cdy = new FXY(ac[3],0.0,ac[2]);
			
				fbx[k] = Utils.interpolateFunctionOnElement(bdx, e);
				fby[k] = Utils.interpolateFunctionOnElement(bdy, e);
				fcx[k] = Utils.interpolateFunctionOnElement(cdx, e);
				fcy[k] = Utils.interpolateFunctionOnElement(cdy, e);
				
				fbxMfcx[k] = fbx[k].M(fcx[k]);
				fbyMfcy[k] = fby[k].M(fcy[k]);
				fbMfc[k] = fb[k].M(fc[k]);
			}
			
			Function integrand = FMath.sum(
						ff,
						fk1.M(FMath.sum(fbxMfcx)),
						fk1.M(FMath.sum(fbyMfcy)),
						fk2.M(FMath.sum(fbMfc))
					).M(v);
			
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
