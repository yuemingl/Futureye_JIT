package edu.uta.futureye.lib.weakform;

import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Face;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.function.FMath;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.function.intf.VectorMathFunc;
import edu.uta.futureye.lib.shapefun.SFConstant0;
import edu.uta.futureye.util.MathEx;
import edu.uta.futureye.util.Utils;
import static edu.uta.futureye.function.FMath.*;

/**
 * <blockquote><pre>
 * Problem:
 * -\nabla{k*\nabla{\vec{u}} 
 * 		+ \vec{U}\cdot\nabla\vec{u}
 * 		+ c*\vec{u} 
 * 		+ \nabla{p} 
 * 		= \vec{f}
 * div{\vec{u}} = 0
 * 
 * Weak form:
 *   find \vec{u} \in H_0^1(div;\Omega), p \in L_2(\Omega)
 *   such that, for all \vec{v} \in H_0^1(div;\Omega), q \in L_2(\Omega)
 *   
 *   (\nabla\vec{v},k*\nabla\vec{u})
 *   		+ (\vec{U}\cdot\nabla\vec{u},\vec{v})
 *   		+ (c*\vec{u},\vec{v})
 *   		- (div{\vec{v}},p)
 *  		+ (q,div{\vec{u}}) 
 *			= (\vec{v},\vec{f})
 *
 *   k* [(v1_x,u1_x) + (v1_y,u1_y) + (v1_z,u1_z) + 
 *       (v2_x,u2_x) + (v2_y,u2_y) + (v2_z,u2_z) +
 *       (v3_x,u3_x) + (v3_y,u3_y) + (v3_z,u3_z) ] +
 *       
 *   	[(U1*u1_x,v1)+(U2*u1_y,v1)+(U3*u1_z,v1)] + 
 *      [(U1*u2_x,v2)+(U2*u2_y,v2)+(U3*u2_z,v2)] +
 *      [(U1*u2_x,v2)+(U2*u2_y,v2)+(U3*u3_z,v3)] +
 *      
 *   		+ c*[(u1,v1)+(u2,v2)+(u3,v3)]
 *   
 *			- (v1_x+v2_y+v3_z,p)
 *			+ (q,u1_x+u2_y+u3_z)
 *
 *			= (v1,f1)+(v2,f2)+(v3,f3)
 *
 * where
 *   \vec{u}=(u1,u2,u3): velocity vector field
 *   \vec{f}=(f1,f2,f3): body force
 *   \vec{U}=(U1,U2,U3): previous velocity
 * </blockquote></pre>
 * 
 * @author liuyueming
 *
 */
public class WeakFormNavierStokes3D extends AbstractVectorWeakForm {
	protected VectorMathFunc g_f = null;
	protected MathFunc g_k = null;
	protected VectorMathFunc g_U = null;
	protected MathFunc g_c = null;
	//Robin:  k*u_n + d*u - p\vec{n} = 0
	protected VectorMathFunc g_d = null;

	public void setF(VectorMathFunc f) {
		this.g_f = f;
	}
	
	public void setParam(MathFunc k, VectorMathFunc U, MathFunc c) {
		this.g_k = k;
		this.g_U = U;
		this.g_c = c;
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
			MathFunc fc = Utils.interpolateOnElement(g_c,e);
			VectorMathFunc fU = Utils.interpolateOnElement(g_U, e);
			
			MathFunc u1 = u.get(1), u2 = u.get(2), u3 = u.get(3), p  = u.get(4);
			MathFunc v1 = v.get(1), v2 = v.get(2), v3 = v.get(3), q  = v.get(4);
			//upwind: v1 and v2 in (u,v,p)
			if(this.testDOF.getVVFComponent() != 4) {
			//if(this.vDOFLocalIndex<=24) {
				Node node1 = testDOF.getNodeOwner();
				//\tidle{v} = v + \tidle{k}*\hat{U}*v,j/\norm{U}
				Vector valU = g_U.value(new Variable().setIndex(node1.globalIndex));
				double normU = valU.norm2();
				Vector valU_hat = valU.scale(1.0/normU);
				double h = e.getElementDiameter();
				double k = g_k.apply(Variable.createFrom(g_k, node1, node1.globalIndex));
				double alpha = normU*h/(2*k);
				//double k_tidle = 2*(MathEx.coth(alpha)-1.0/alpha)*normU*h;
				double k_tidle = (MathEx.coth(alpha)-1.0/alpha)*normU*h;
				if(!(v1.isConstant())) {
					v1.A(FMath.grad(v1,"x","y","z").dot(valU_hat).M(k_tidle).D(valU.norm2()));
				}
				if(!(v2.isConstant())) {
					v2.A(FMath.grad(v2,"x","y","z").dot(valU_hat).M(k_tidle).D(valU.norm2()));
				}
			}
			
			
/**
 *   k* [(v1_x,u1_x) + (v1_y,u1_y) + (v1_z,u1_z) + 
 *       (v2_x,u2_x) + (v2_y,u2_y) + (v2_z,u2_z) +
 *       (v3_x,u3_x) + (v3_y,u3_y) + (v3_z,u3_z) ] +
 *       
 *   	[(U1*u1_x,v1)+(U2*u1_y,v1)+(U3*u1_z,v1)] + 
 *      [(U1*u2_x,v2)+(U2*u2_y,v2)+(U3*u2_z,v2)] +
 *      [(U1*u2_x,v2)+(U2*u2_y,v2)+(U3*u3_z,v3)] +
 *      
 *   		+ c*[(u1,v1)+(u2,v2)+(u3,v3)]
 *   
 *			- (v1_x+v2_y+v3_z,p)
 *			+ (q,u1_x+u2_y+u3_z)
 *
 *			= (v1,f1)+(v2,f2)+(v3,f3)
 * 
 */
			VectorMathFunc grad_u1 = grad(u1,"x","y","z");
			VectorMathFunc grad_u2 = grad(u2,"x","y","z");
			VectorMathFunc grad_u3 = grad(u3,"x","y","z");
			MathFunc uv1 = grad_u1.dot( grad(v1,"x","y","z") );
			MathFunc uv2 = grad_u2.dot( grad(v2,"x","y","z") );
			MathFunc uv3 = grad_u3.dot( grad(v3,"x","y","z") );
			MathFunc div_v = v1.diff("x").A(v2.diff("y")).A(v3.diff("z"));
			MathFunc div_u = u1.diff("x").A(u2.diff("y")).A(u3.diff("z"));
			MathFunc cvect = fU.dot(grad_u1).M(v1).A(fU.dot(grad_u2).M(v2)).A(fU.dot(grad_u3).M(v3));
			MathFunc cuv = fc.M(u1.M(v1).A(u2.M(v2)).A(u3.M(v3)));
			integrand = fk.M( uv1.A(uv2).A(uv3) ).A( cvect ).A( cuv ).S( div_v.M(p) ).A( div_u.M(q) );
			return integrand;
		} else if(itemType==ItemType.Border) {
			if(g_d != null) {
				Element be = e;
				MathFunc fd1 = Utils.interpolateOnElement(g_d.get(1), be);
				MathFunc fd2 = Utils.interpolateOnElement(g_d.get(2), be);
				MathFunc fd3 = Utils.interpolateOnElement(g_d.get(3), be);
				MathFunc u1 = u.get(1), u2 = u.get(2), u3 = u.get(3);
				MathFunc v1 = v.get(1), v2 = v.get(2), v3 = v.get(3);
				//Robin:  - k*u_n = d*u - p\vec{n}
				//d1*u1 + d2*u2 + d3*u3
				MathFunc borderIntegrand = fd1.M(u1.M(v1)).A(fd2.M(u2.M(v2))).A(fd3.M(u3.M(v3)));
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
			MathFunc f3 = Utils.interpolateOnElement(g_f.get(3), e);
			MathFunc v1 = v.get(1), v2 = v.get(2), v3 = v.get(3);
			//(v1*f1+v2*f2)
			MathFunc integrand = v1.M(f1).A(v2.M(f2)).A(v3.M(f3));
			return integrand;
		} else if(itemType==ItemType.Border) {
			Element be = e;
			MathFunc v1 = v.get(1), v2 = v.get(2), v3 = v.get(3), p  = v.get(4);
			//Robin:  - k*u_n = d*u - p\vec{n}
			//- p\vec{n} = - p*n1*v1 - p*n2*v2 - p*n3*v3
			Face face = (Face)be.getGeoEntity();
			Vector n = face.getNormVector();
			MathFunc n1 = C(-1.0*n.get(1)), n2 = C(-1.0*n.get(2)), n3 = C(-1.0*n.get(3));
			MathFunc borderIntegrand = p.M(v1.M(n1)).A(p.M(v2.M(n2))).A(p.M(v3.M(n3)));
			return borderIntegrand;
		}
		return null;
	}
	
	public boolean isVVFComponentCoupled(int nComponent1, int nComponent2) {
		if(nComponent1 == nComponent2) return true;
		else if(nComponent1 == 4 || nComponent2 == 4) return true;
		else return false;
	}
}
