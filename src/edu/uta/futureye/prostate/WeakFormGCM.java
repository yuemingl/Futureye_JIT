package edu.uta.futureye.prostate;

import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.intf.WeakForm;
import edu.uta.futureye.function.DerivativeIndicator;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.function.intf.ShapeFunction;
import edu.uta.futureye.function.operator.FOBasic;
import edu.uta.futureye.function.operator.FOIntegrate;
import edu.uta.futureye.util.Utils;

/**
 * Solve
 *  (k*u_x, v_x) + (k*u_y, v_y) + (b1*u_x, v) + (b2*u_y, v) + (c*u, v)= (f, v)
 * 
 * @author liuyueming
 *
 */
public class WeakFormGCM implements WeakForm {

	protected ShapeFunction u = null;
	protected ShapeFunction v = null;
	protected int uDOFLocalIndex;
	protected int vDOFLocalIndex;
	
	protected Function g_f = null;
	
	protected Function g_k = null;
	protected Function g_c = null;
	protected Function g_b1 = null;
	protected Function g_b2 = null;
	
	protected Function g_q = null;
	protected Function g_d = null;

	public void setF(Function f) {
		this.g_f = f;
	}
	
	public void setParam(Function k,Function c,Function b1,Function b2) {
		this.g_k = k;
		this.g_c = c;
		this.g_b1 = b1;
		this.g_b2 = b2;
	}
	
	//Robin: d*u + k*u_n = q
	public void setRobin(Function q,Function d) {
		this.g_q = q;
		this.g_d = d;
	}	

	@Override
	public Function leftHandSide(Element e, ItemType itemType) {
		u.asignElement(e);
		v.asignElement(e);

		if(itemType==ItemType.Domain)  {
			//Integrand part of Weak Form on element e
			DerivativeIndicator di_x = new DerivativeIndicator("x1");
			DerivativeIndicator di_y = new DerivativeIndicator("y1");
			Function integrand = null;
				
			Function fk = Utils.interplateFunctionOnElement(g_k,e);
			Function fc = Utils.interplateFunctionOnElement(g_c,e);
			Function fb1 = Utils.interplateFunctionOnElement(g_b1,e);
			Function fb2 = Utils.interplateFunctionOnElement(g_b2,e);
			
			integrand = FOBasic.PlusAll(
						FOBasic.Mult(fk, FOBasic.Plus(
								FOBasic.Mult(u.derivative(di_x), v.derivative(di_x)),
								FOBasic.Mult(u.derivative(di_y), v.derivative(di_y))
						)),
						FOBasic.Mult(fb1,FOBasic.Mult(u.derivative(di_x), v)),
						FOBasic.Mult(fb2,FOBasic.Mult(u.derivative(di_y), v)),
						FOBasic.Mult(fc, FOBasic.Mult(u, v))
					);

			//Numerical integration on element e
			Function integral = null;
			if(e.getVertexList().size() == 3) {
				integral = FOIntegrate.intOnTriangleRefElement(
						FOBasic.Mult(integrand, e.getJacobin()),5
						);
			} else if (e.getVertexList().size() == 4) {
				integral = FOIntegrate.intOnRectangleRefElement(
						FOBasic.Mult(integrand, e.getJacobin()),2 //TODO !!! 1
						);
			}
	
			return integral;
		}
		else if(itemType==ItemType.Border) {
			if(g_d != null) {
				Element be = e;
				Function fd = Utils.interplateFunctionOnElement(g_d, e);
				Function borderIntegrand = FOBasic.Mult(FOBasic.Mult(fd, u), v);
				//Numerical integration the border of element e
				//System.out.println(be.getJacobin().value(null));
				Function borderIntgral = FOIntegrate.intOnLinearRefElement(
						FOBasic.Mult(borderIntegrand, be.getJacobin()),5
					);
				return borderIntgral;
			}
		}
		return null;
	}

	@Override
	public Function rightHandSide(Element e, ItemType itemType) {
		v.asignElement(e);
		if(itemType==ItemType.Domain)  {
			Function ff = Utils.interplateFunctionOnElement(g_f, e);
			Function integrand = FOBasic.Mult(ff,v);
			
			Function integral = null;
			if(e.getVertexList().size() == 3) {
				integral = FOIntegrate.intOnTriangleRefElement(
						FOBasic.Mult(integrand, e.getJacobin()),5
					);
			} else if (e.getVertexList().size() == 4) {
				integral = FOIntegrate.intOnRectangleRefElement(
						FOBasic.Mult(integrand, e.getJacobin()),2 //TODO !!! 1
						);
			}
			//System.out.println(e.nodes.size()+":  "+integral.value(null));
			return integral;
		} else if(itemType==ItemType.Border) {//add 2010-9-28
			Element be = e;
			Function fq = Utils.interplateFunctionOnElement(g_q, e);
			Function borderIntegrand = FOBasic.Mult(fq, v);
			//Numerical integration the border of element e
			//System.out.println(be.getJacobin().value(null));
		
			Function borderIntgral = FOIntegrate.intOnLinearRefElement(
					FOBasic.Mult(borderIntegrand, be.getJacobin()),5
				);
			//System.out.println(borderIntgral.value(null));
			return borderIntgral;
		} else
			return null;	
	}

	/**
	 * N = N(r,s) = N( r(x,y), s(x,y) )
	 */
	@Override
	public void setShapeFunction(ShapeFunction trial, int trialDofLocalIndex,
			ShapeFunction test, int testDofLocalIndex) {
		u = trial;
		v = test;
		this.uDOFLocalIndex = trialDofLocalIndex;
		this.vDOFLocalIndex = testDofLocalIndex;
	}


}
