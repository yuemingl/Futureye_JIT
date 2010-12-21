package edu.uta.futureye.core;

import java.util.List;

import edu.uta.futureye.algebra.Matrix;
import edu.uta.futureye.core.intf.WeakForm;
import edu.uta.futureye.function.DerivativeIndicator;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.FConstant;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.function.intf.FunctionDerivable;
import edu.uta.futureye.function.intf.ShapeFunction;
import edu.uta.futureye.function.operator.FOBasic;
import edu.uta.futureye.function.operator.FOIntegrate;
import edu.uta.futureye.test.VectorBasedFunction;
import edu.uta.futureye.util.DOFList;
import edu.uta.futureye.util.PairElementMatrix;
import edu.uta.futureye.util.Utils;

/**
 * 系数反问题
 * Solve: (U*u, v) = (f, v) - (k*grad(U),grad(v))
 * where u is unknown
 * U,f and k is known
 * 
 * @author liuyueming
 */
public class WeakFormL22D implements WeakForm {
	protected ShapeFunction u = null;
	protected ShapeFunction v = null;
	protected int uDOFLocalIndex;
	protected int vDOFLocalIndex;
	protected VectorBasedFunction g_U = null;
	protected Function g_f = null;
	protected Function g_k = null;


	@Override
	public Function leftHandSide(Element e, ItemType itemType) {
		u.asignElement(e);
		v.asignElement(e);

		if(itemType==ItemType.Domain)  {
			//Integrand part of Weak Form on element e
			Function integrand = null;
	
//			int nNode = e.nodes.size();
//			FunctionDerivable fU = new FConstant(0.0);
//			for(int j=1;j<=nNode;j++) {
//				DOFList dofListI = e.getDOFList(j);
//				for(int k=1;k<=dofListI.size();k++) {
//					DOF dofI = dofListI.at(k);
//					FunctionDerivable UValue = new FConstant(g_U.value(new Variable(dofI.globalIndex)));
//						fU = FOBasicDerivable.Plus(fU, 
//								FOBasicDerivable.Mult(UValue, dofI.shapeFunction));
//				}
//			}
			
			Function fU = Utils.interplateFunctionOnElement(g_U, e);
			integrand = FOBasic.Mult(fU, FOBasic.Mult(u, v));
	
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
		return null;
	}

	@Override
	public Function rightHandSide(Element e, ItemType itemType) {
		v.asignElement(e);
		if(itemType==ItemType.Domain)  {
			Function ff = Utils.interplateFunctionOnElement(g_f, e);
			Function fk = Utils.interplateFunctionOnElement(g_k, e);
			
			DerivativeIndicator di_x = new DerivativeIndicator("x1");
			DerivativeIndicator di_y = new DerivativeIndicator("y1");
			Function rlt_dx = new FConstant(0.0);
			Function rlt_dy = new FConstant(0.0);
			int nNode = e.nodes.size();
			for(int i=1;i<=nNode;i++) {
				DOFList dofListI = e.getDOFList(i);
				for(int k=1;k<=dofListI.size();k++) {
					DOF dofI = dofListI.at(k);
					Variable var = Variable.createFrom(g_U, dofI.getOwnerNode(), dofI.getGlobalNumber());
					FunctionDerivable PValue = new FConstant(g_U.value(var));
					ShapeFunction shape = dofI.getShapeFunction();
					//TODO 不要忘记！！！
					shape.asignElement(e);
					rlt_dx = FOBasic.Plus(rlt_dx, FOBasic.Mult(PValue, shape.derivative(di_x)));
					rlt_dy = FOBasic.Plus(rlt_dy, FOBasic.Mult(PValue, shape.derivative(di_y)));
				}
			}
			
			Function integrand = FOBasic.Minus(
				FOBasic.Mult(ff,v),
				FOBasic.Mult(fk,
						FOBasic.Plus(
								FOBasic.Mult(rlt_dx, v.derivative(di_x)),
								FOBasic.Mult(rlt_dy, v.derivative(di_y))
			)));
			
			
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
		}
		return null;
	}

	@Override
	public void setShapeFunction(ShapeFunction trial, int trialDofLocalIndex,
			ShapeFunction test, int testDofLocalIndex) {
		u = trial;
		v = test;
		this.uDOFLocalIndex = trialDofLocalIndex;
		this.vDOFLocalIndex = testDofLocalIndex;
	}
	
	public void setF(Function f) {
		this.g_f = f;
	}
	
	public void setParam(Function k,VectorBasedFunction U) {
		this.g_k = k;
		this.g_U = U;
	}

	@Override
	public List<PairElementMatrix> associateElement(Element e) {
		// TODO Auto-generated method stub
		return null;
	}
}
