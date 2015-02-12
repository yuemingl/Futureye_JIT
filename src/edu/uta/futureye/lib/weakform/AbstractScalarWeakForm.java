package edu.uta.futureye.lib.weakform;

import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.core.DOF;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.intf.WeakForm;
import edu.uta.futureye.function.FEMFunc;
import edu.uta.futureye.function.intf.MathFun;
import edu.uta.futureye.function.intf.ScalarShapeFunction;
import edu.uta.futureye.function.operator.FOIntegrate;
import edu.uta.futureye.util.FutureyeException;

public abstract class AbstractScalarWeakForm implements WeakForm {
	protected DOF trialDOF = null;
	protected DOF testDOF = null;
	protected ScalarShapeFunction u = null;
	protected ScalarShapeFunction v = null;
	protected int uDOFLocalIndex;
	protected int vDOFLocalIndex;

	@Override
	public void assembleElement(Element e, 
			Matrix globalStiff,	Vector globalLoad) {
		throw new UnsupportedOperationException();
	}

	@Override
	public FEMFunc leftHandSide(Element e, ItemType itemType) {
		throw new UnsupportedOperationException();
	}

	@Override
	public FEMFunc rightHandSide(Element e, ItemType itemType) {
		throw new UnsupportedOperationException();
	}

//	@Override
//	public void setShapeFunction(ShapeFunction trial, int trialDofLocalIndex,
//			ShapeFunction test, int testDofLocalIndex) {
//		this.u = (ScalarShapeFunction)trial;
//		this.v = (ScalarShapeFunction)test;
//		this.uDOFLocalIndex = trialDofLocalIndex;
//		this.vDOFLocalIndex = testDofLocalIndex;
//	}

	@Override
	public void setDOF(DOF trialDOF, DOF testDOF) {
		this.trialDOF = trialDOF;
		this.testDOF = testDOF;
		if(trialDOF != null) {
			this.u = trialDOF.getSSF();
			this.uDOFLocalIndex = trialDOF.getLocalIndex();
		}
		if(testDOF != null) {
			this.v = testDOF.getSSF();
			this.vDOFLocalIndex = testDOF.getLocalIndex();
		}
	}
	
	@Override
	public DOF getTrialDOF() {
		return this.trialDOF;
	}
	
	@Override
	public DOF getTestDOF() {
		return this.testDOF;
	}
	
	@Override
	public double integrate(Element e, MathFun fun) {
		if(fun == null) return 0.0;
		if(e.dim() == 2) {
			if(e.vertices().size() == 3) {
				//三角形单元
				return FOIntegrate.intOnTriangleRefElement(
						fun.M(e.getJacobin()),2
						);
			} else if (e.vertices().size() == 4) {
				//四边形单元
				return FOIntegrate.intOnRectangleRefElement(
						fun.M(e.getJacobin()),5 //TODO
						);
			}
		} else if(e.dim() == 3) {
			if(e.vertices().size() == 4) {
				//四面体单元
				return FOIntegrate.intOnTetrahedraRefElement(
						fun.M(e.getJacobin()),2
					);
			} else if(e.vertices().size() == 8) {
				//六面体单元
				return FOIntegrate.intOnHexahedraRefElement(
						fun.M(e.getJacobin()),2
						
					);
			}
		} else if(e.dim() == 1) {
			//一维单元
			return FOIntegrate.intOnLinearRefElement(
					fun.M(e.getJacobin()),5
				);
		} else {
			throw new FutureyeException(
					"Can NOT integrate on e" + e.vertices());
		}
		throw new FutureyeException("Error");
	}
	
	/**
	 * No meaning for scalar valued problems
	 */
	public boolean isVVFComponentCoupled(int nComponent1, int nComponent2) {
		return false;
	}
	
	public void preProcess(Element e) {
	}
}
