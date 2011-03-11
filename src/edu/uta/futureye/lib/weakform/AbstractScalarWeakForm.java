package edu.uta.futureye.lib.weakform;

import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.intf.WeakForm;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.function.intf.ScalarShapeFunction;
import edu.uta.futureye.function.intf.ShapeFunction;
import edu.uta.futureye.function.operator.FOBasic;
import edu.uta.futureye.function.operator.FOIntegrate;
import edu.uta.futureye.util.FutureyeException;

public class AbstractScalarWeakForm implements WeakForm {
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
	public Function leftHandSide(Element e, ItemType itemType) {
		throw new UnsupportedOperationException();
	}

	@Override
	public Function rightHandSide(Element e, ItemType itemType) {
		throw new UnsupportedOperationException();
	}

	@Override
	public void setShapeFunction(ShapeFunction trial, int trialDofLocalIndex,
			ShapeFunction test, int testDofLocalIndex) {
		this.u = (ScalarShapeFunction)trial;
		this.v = (ScalarShapeFunction)test;
		this.uDOFLocalIndex = trialDofLocalIndex;
		this.vDOFLocalIndex = testDofLocalIndex;
	}

	@Override
	public double integrate(Element e, Function fun) {
		if(fun == null) return 0.0;
		Function integral = null;
		if(e.eleDim() == 2) {
			if(e.vertices().size() == 3) {
				//三角形单元
				integral = FOIntegrate.intOnTriangleRefElement(
						FOBasic.Mult(fun, e.getJacobin()),5
						);
			} else if (e.vertices().size() == 4) {
				//四边形单元
				integral = FOIntegrate.intOnRectangleRefElement(
						FOBasic.Mult(fun, e.getJacobin()),2 //TODO
						);
			}
		} else if(e.eleDim() == 3) {
			if(e.vertices().size() == 4) {
				//四面体单元
				integral = FOIntegrate.intOnTetrahedraRefElement(
						FOBasic.Mult(fun, e.getJacobin()),2
					);
			}
		} else if(e.eleDim() == 1) {
			//一维单元
			integral = FOIntegrate.intOnLinearRefElement(
					FOBasic.Mult(fun, e.getJacobin()),5
				);
		} else {
			FutureyeException ex = new FutureyeException(
					"Can integrate on e" + e.vertices());
			ex.printStackTrace();
			System.exit(-1);
		}
		return integral.value(null);
	}
}
