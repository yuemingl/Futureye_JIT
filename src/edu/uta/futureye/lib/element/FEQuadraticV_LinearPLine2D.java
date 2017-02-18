package edu.uta.futureye.lib.element;

import static edu.uta.futureye.function.FMath.C0;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Line2DCoord;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.intf.CoordTrans;
import edu.uta.futureye.core.intf.VecFiniteElement;
import edu.uta.futureye.function.basic.FX;
import edu.uta.futureye.function.basic.SpaceVectorFunction;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.function.intf.VecMathFunc;

/**
 * P2/P1
 * -Continuous quadratic velocity
 * -Piecewise linear pressure
 * 
 * Restrict to boundary edge:
 * * Velocity:
 * | \
 * |  \
 * |   \
 * |    \
 * |     \
 * 1--3-- 2
 * NV1,NV2,NV3
 * 
 * NV1 = r*(r-1)/2
 * NV2 = (r+1)*r/2
 * NV3 = 1-r*r
 * 
 * * Pressure:
 * | \
 * |  \
 * |   \
 * |    \
 * 1---- 2 
 * NP1,NP2
 * 
 * NP1 = (1-r)/2
 * NP2 = (1+r)/2
 * 
 * Vector valued shape functions
 * Ni = (v1,v2,p)', i=1,...,8, on boundary
 * where
 * N1 =  (NV1,0,0)'
 * N2 =  (NV2,0,0)'
 * N3 =  (NV3,0,0)'
 * N4 =  (0,NV1,0)'
 * N5 =  (0,NV2,0)'
 * N6 =  (0,NV3,0)'
 * N7 =  (0,0,NP1)'
 * N8 =  (0,0,NP2)'
 * 
 * @param funIndex
 * @return
 */


public class FEQuadraticV_LinearPLine2D implements VecFiniteElement {
	Line2DCoord coord;
	
	//Construct a function with the coordinate of points in an element as parameters
	String[] argsOrder;
	
	public int nDOFs = 2+2+1;
	VecMathFunc[] shapeFuncs = new VecMathFunc[nDOFs];

	public FEQuadraticV_LinearPLine2D() {
		FX x1 = new FX("x1");
		FX x2 = new FX("x2");
		FX y1 = new FX("y1");
		FX y2 = new FX("y2");
		
		MathFunc r = coord.getCoordR();
		
		this.argsOrder = new String[]{x1,x2,y1,y2,r};
		
		
		//shape function for velocity component
		MathFunc NV1 = r*(r-1)/2;
		MathFunc NV2 = (r+1)*r/2;
		MathFunc NV3 = 1-r*r;
		
		//shape functions for pressure component
		MathFunc NP1 = (1-r)/2;
		MathFunc NP2 = (1+r)/2;

		shapeFuncs[0] = new SpaceVectorFunction(NV1, C0, C0);
		shapeFuncs[1] = new SpaceVectorFunction(NV2, C0, C0);
		shapeFuncs[2] = new SpaceVectorFunction(NV3, C0, C0);
		shapeFuncs[3] = new SpaceVectorFunction(C0, NV1, C0);
		shapeFuncs[4] = new SpaceVectorFunction(C0, NV2, C0);
		shapeFuncs[5] = new SpaceVectorFunction(C0, NV3, C0);
		shapeFuncs[6] = new SpaceVectorFunction(C0, C0, NP1);
		shapeFuncs[7] = new SpaceVectorFunction(C0, C0, NP2);
	}

	@Override
	public int getNumberOfDOFs() {
		return this.nDOFs;
	}

	@Override
	public VecMathFunc[] getShapeFunctions() {
		return this.shapeFuncs;
	}

	@Override
	public String[] getArgsOrder() {
		return this.argsOrder;
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
	
	@Override
	public int getGlobalIndex(Mesh mesh, Element e, int localIndex) {
		if(localIndex>=1 && localIndex <= 3) {
			return e.nodes.at(localIndex).globalIndex;
		} else if(localIndex>=4 && localIndex<=6) {
			int nNode = mesh.getNodeList().size();
			return nNode + e.nodes.at(localIndex-3).globalIndex;
		} else if(localIndex>=7 && localIndex<=8) {
			int nNode = mesh.getNodeList().size();
			return 2*nNode + e.vertices().at(localIndex-8).globalNode().globalIndex;
		} else {
			throw new RuntimeException("local index = "+localIndex+". It should be in 1...8");
		}
	}

	@Override
	public int getTotalNumberOfDOFs(Mesh mesh) {
		throw new UnsupportedOperationException("Call FEQuadraticV_LinearP.getTotalNumberOfDOFs() intstead");
	}

	@Override
	public int getVVFComponentIndex(int localIndex) {
		if(localIndex >= 1 && localIndex <= 3)
			return 1;
		else if(localIndex >= 4 && localIndex <= 6)
			return 2;
		else if(localIndex ==7 || localIndex == 8)
			return 3;
		else
			throw new RuntimeException("local index should be in the range of [1,"+(shapeFuncs.length+1)+"]");
	}
	
	@Override
	public CoordTrans getCoordTrans() {
		return this.coord;
	}
}
