package edu.uta.futureye.lib.weakform;

import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.core.Edge;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.function.intf.*;
import edu.uta.futureye.util.Utils;
import static edu.uta.futureye.function.FMath.*;

/**
 * Problem:
 * -\nabla{k*\nabla{\vec{u}}} + \nabla{p} = \vec{f}
 * div{\vec{u}} = 0
 * 
 * Weak form:
 *   find \vec{u} \in H_0^1(div;\Omega), p \in L_2(\Omega)
 *   such that, for all \vec{v} \in H_0^1(div;\Omega), q \in L_2(\Omega)
 *   
 *   (\nabla{\vec{v}},k*\nabla{\vec{u}}) - (div{\vec{v}},p) 
 *                   + (q,div{\vec{u}}) = (\vec{v},\vec{f})
 *
 *   (v1_x,k*u1_x) + (v1_y,k*u1_y) + (v2_x,k*u2_x) + (v2_y,k*u2_y) 
 *                   - (v1_x+v2_y,p) + (q,u1_x+u2_y) = (v1,f1)+(v2,f2)    
 *
 * where
 *   \vec{u}=(u1,u2): velocity vector field    
 *   \vec{f}=(f1,f2): body force
 *   
 * @author liuyueming
 *
 */
public class WeakFormStokes extends AbstractVectorWeakForm {
	protected VectorMathFunc g_f = null;
	protected MathFunc g_k = null;
	//Robin:  k*u_n + d*u - p\vec{n} = 0
	protected VectorMathFunc g_d = null;

	public void setF(VectorMathFunc f) {
		this.g_f = f;
	}
	
	public void setParam(MathFunc k) {
		this.g_k = k;
	}
	
	//Robin:  k*u_n + d*u - p\vec{n} = 0
	public void setRobin(VectorMathFunc d) {
		this.g_d = d;
	}
	
	@Override
	public MathFunc leftHandSide(Element e, ItemType itemType) {

		if(itemType==ItemType.Domain)  {
			//Integrand part of Weak Form on element e
			MathFunc integrand = null;
			MathFunc fk = Utils.interpolateOnElement(g_k,e);
			MathFunc u1 = u.get(1), u2 = u.get(2), p  = u.get(3);
			MathFunc v1 = v.get(1), v2 = v.get(2), q  = v.get(3);
			//(v1_x,k*u1_x) + (v1_y,k*u1_y) + 
			//(v2_x,k*u2_x) + (v2_y,k*u2_y) - 
			//(v1_x+v2_y,p) + (q,u1_x+u2_y) 
			MathFunc uv1 = grad(u1,"x","y" ).dot( grad(v1,"x","y") );
			MathFunc uv2 = grad(u2,"x","y" ).dot( grad(v2,"x","y") );
			MathFunc div_v = v1.diff("x").A(v2.diff("y"));
			MathFunc div_u = u1.diff("x").A(u2.diff("y"));
			integrand = fk.M( uv1.A(uv2) ).S( div_v.M(p) ).A( div_u.M(q) );
			return integrand;
		}
		else if(itemType==ItemType.Border) {
			if(g_d != null) {
				Element be = e;
				MathFunc fd1 = Utils.interpolateOnElement(g_d.get(1), be);
				MathFunc fd2 = Utils.interpolateOnElement(g_d.get(2), be);
				MathFunc u1 = u.get(1), u2 = u.get(2);
				MathFunc v1 = v.get(1), v2 = v.get(2);
				//Robin:  - k*u_n = d*u - p\vec{n}
				//d1*u1 + d2*u2
				MathFunc borderIntegrand = fd1.M(u1.M(v1)).A(fd2.M(u2.M(v2)));
				return borderIntegrand;
			}
		}
		return null;
	}

	@Override
	public MathFunc rightHandSide(Element e, ItemType itemType) {
		if(itemType==ItemType.Domain)  {
			MathFunc f1 = Utils.interpolateOnElement(g_f.get(1), e);
			MathFunc f2 = Utils.interpolateOnElement(g_f.get(2), e);
			MathFunc v1 = v.get(1);
			MathFunc v2 = v.get(2);
			//(v1*f1+v2*f2)
			MathFunc integrand = v1.M(f1).A(v2.M(f2));
			return integrand;
		} else if(itemType==ItemType.Border) {
			Element be = e;
			MathFunc v1 = v.get(1), v2 = v.get(2), p  = v.get(3);
			//Robin:  - k*u_n = d*u - p\vec{n}
			//- p\vec{n} = - p*n1*v1 - p*n2*v2
			Edge edge = (Edge)be.getGeoEntity();
			Vector n = edge.getNormVector();
			MathFunc n1 = C(-1.0*n.get(1));
			MathFunc n2 = C(-1.0*n.get(2));
			MathFunc borderIntegrand = p.M(v1.M(n1)).A(p.M(v2.M(n2)));
			return borderIntegrand;
		}
		return null;
	}
	
	public boolean isVVFComponentCoupled(int nComponent1, int nComponent2) {
		if(nComponent1 == nComponent2) return true;
		else if(nComponent1 == 3 || nComponent2 == 3) return true;
		else return false;
	}
}
