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

/**
 * 用有限元方法求导数
 * Solve: (w, v) = (U_x, v)
 * where w is unknown
 *   U_x is the piecewise derivative on mesh
 *   w is an approximation of U_x
 *   
 * @author liuyueming
 */
public class WeakFormDerivative implements WeakForm {
	protected ShapeFunction u = null;
	protected ShapeFunction v = null;
	protected int uDOFLocalIndex;
	protected int vDOFLocalIndex;
	
	protected VectorBasedFunction g_U = null;
	protected String varName; // "x" or "y"

	public WeakFormDerivative(String varName) {
		this.varName = varName;
	}
	
	@Override
	public Function leftHandSide(Element e, ItemType itemType) {
		u.asignElement(e);
		v.asignElement(e);

		if(itemType==ItemType.Domain)  {
			//Integrand part of Weak Form on element e
			Function integrand = null;
	
			integrand = FOBasic.Mult(u, v);
	
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
			
			DerivativeIndicator di = new DerivativeIndicator(varName+"1");
			Function rlt = new FConstant(0.0);
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
					rlt = FOBasic.Plus(rlt, FOBasic.Mult(PValue, shape.derivative(di)));
				}
			}
			
			Function integrand = FOBasic.Mult(rlt, v);
			
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
	
	public void setParam(VectorBasedFunction U) {
		this.g_U = U;
	}

	@Override
	public List<PairElementMatrix> associateElement(Element e) {
		// TODO Auto-generated method stub
		return null;
	}
}
