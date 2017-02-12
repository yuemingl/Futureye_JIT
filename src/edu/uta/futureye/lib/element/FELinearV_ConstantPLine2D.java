package edu.uta.futureye.lib.element;

import static edu.uta.futureye.function.FMath.C0;
import static edu.uta.futureye.function.FMath.C1;

import java.util.Map;

import edu.uta.futureye.core.Line2DCoord;
import edu.uta.futureye.core.intf.VecFiniteElement;
import edu.uta.futureye.function.basic.FX;
import edu.uta.futureye.function.basic.SpaceVectorFunction;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.function.intf.VectorMathFunc;

public class FELinearV_ConstantPLine2D implements VecFiniteElement {
	Line2DCoord coord;
	
	//Construct a function with the coordinate of points in an element as parameters
	String[] argsOrder;
	
	public int nDOFs = 2+2+1;
	VectorMathFunc[] shapeFuncs = new VectorMathFunc[nDOFs];

	public FELinearV_ConstantPLine2D() {
		FX x1 = new FX("x1");
		FX x2 = new FX("x2");
		FX y1 = new FX("y1");
		FX y2 = new FX("y2");

		this.coord = new Line2DCoord(x1,x2,y1,y2);
		
		MathFunc r = this.coord.getCoordR();
		this.argsOrder = new String[]{x1,x2,y1,y2,r};
		
		//shape function for velocity component
		MathFunc NV1 = 0.5*(1-r);
		MathFunc NV2 = 0.5*(1+r);

		shapeFuncs[0] = new SpaceVectorFunction(NV1, C0, C0);
		shapeFuncs[1] = new SpaceVectorFunction(NV2, C0, C0);
		shapeFuncs[2] = new SpaceVectorFunction(C0, NV1, C0);
		shapeFuncs[3] = new SpaceVectorFunction(C0, NV2, C0);
		shapeFuncs[4] = new SpaceVectorFunction(C0, C0, C1);
	}

	@Override
	public int getNumberOfDOFs() {
		return this.nDOFs;
	}

	@Override
	public VectorMathFunc[] getShapeFunctions() {
		return this.shapeFuncs;
	}

	@Override
	public Map<String, MathFunc> getCoordTransMap() {
		return this.coord.getCoordTransMap();
	}

	@Override
	public String[] getArgsOrder() {
		return this.argsOrder;
	}

	@Override
	public MathFunc getJacobian() {
		return this.coord.getJacobian();
	}

	@Override
	public VecFiniteElement getBoundaryFE() {
		return null;
	}

	@Override
	public boolean isDOFCoupled(int idx1, int idx2) {
		if(idx1 == 4 || idx2 == 4)
			return true;
		if(idx1 <= 1 && idx2 >=2)
			return false;
		if(idx2 <= 1 && idx1 >=2)
			return false;
		return true;
	}

}
