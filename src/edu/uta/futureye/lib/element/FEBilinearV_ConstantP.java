package edu.uta.futureye.lib.element;


import static edu.uta.futureye.function.FMath.C0;
import static edu.uta.futureye.function.FMath.C1;

import java.util.Map;

import edu.uta.futureye.core.RectAreaCoord;
import edu.uta.futureye.core.intf.VecFiniteElement;
import edu.uta.futureye.function.basic.FX;
import edu.uta.futureye.function.basic.SpaceVectorFunction;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.function.intf.VectorMathFunc;

/**
 * 2D Q1/P0 Element
 * -Continuous bilinear velocity
 * -Piecewise constant pressure
 * 
 * * Velocity: Bilinear shape function: SFBilinearLocal2D
 * * 速度：四边形局部坐标，双线性函数
 * 
 * 4----3
 * |    |
 * |    |
 * 1----2
 * 
 * NV = NV(r,s) = NV( r(x,y), s(x,y) )
 * NV1 = (1-r)*(1-s)/4
 * NV2 = (1+r)*(1-s)/4
 * NV3 = (1+r)*(1+s)/4
 * NV4 = (1-r)*(1+s)/4
 * 
 * * Pressure: Piecewise constant shape function: SFConstant1
 * * 压强：分片常数型函数
 * NP=1
 * 
 * * 2D vector valued shape functions
 * * 二维单元上的形函数，速度压强共9个自由度：
 * Ni = (u1,u2,p)', i=1,...,9
 * 
 * N1  =  (NV1, 0, 0)'
 * N2  =  (NV2, 0, 0)'
 * N3  =  (NV3, 0, 0)'
 * N4  =  (NV4, 0, 0)'
 * N5  =  (0, NV1, 0)'
 * N6  =  (0, NV2, 0)'
 * N7  =  (0, NV3, 0)'
 * N8 =   (0, NV4, 0)'
 * N9 =   (0, 0, NP)'
 *
 */
public class FEBilinearV_ConstantP implements VecFiniteElement {
	RectAreaCoord coord;
	
	//Construct a function with the coordinate of points in an element as parameters
	String[] argsOrder;
	
	public int nDOFs = 4+4+1;
	VectorMathFunc[] shapeFuncs = new VectorMathFunc[nDOFs];

	public FEBilinearV_ConstantP() {
		FX x1 = new FX("x1");
		FX x2 = new FX("x2");
		FX x3 = new FX("x3");
		FX x4 = new FX("x4");
		FX y1 = new FX("y1");
		FX y2 = new FX("y2");
		FX y3 = new FX("y3");
		FX y4 = new FX("y4");
		
		this.coord = new RectAreaCoord(x1,x2,x3,x4,y1,y2,y3,y4);
		
		MathFunc r = coord.getCoordR();
		MathFunc s = coord.getCoordS();
		
		this.argsOrder = new String[]{x1,x2,x3,x4,y1,y2,y3,y4,r,s};
		
		MathFunc NV1 = (1-r)*(1-s)/4;
		MathFunc NV2 = (1+r)*(1-s)/4;
		MathFunc NV3 = (1+r)*(1+s)/4;
		MathFunc NV4 = (1-r)*(1+s)/4;

		shapeFuncs[0] = new SpaceVectorFunction(NV1, C0, C0);
		shapeFuncs[1] = new SpaceVectorFunction(NV2, C0, C0);
		shapeFuncs[2] = new SpaceVectorFunction(NV3, C0, C0);
		shapeFuncs[3] = new SpaceVectorFunction(NV4, C0, C0);
		shapeFuncs[4] = new SpaceVectorFunction(C0, NV1, C0);
		shapeFuncs[5] = new SpaceVectorFunction(C0, NV2, C0);
		shapeFuncs[6] = new SpaceVectorFunction(C0, NV3, C0);
		shapeFuncs[7] = new SpaceVectorFunction(C0, NV4, C0);
		shapeFuncs[8] = new SpaceVectorFunction(C0, C0, C1);
	}

	@Override
	public VectorMathFunc[] getShapeFunctions() {
		return this.shapeFuncs;
	}

	@Override
	public int getNumberOfDOFs() {
		return this.nDOFs;
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

//	//???
//	// we need total number of nodes in the mesh to assign global index for u[2]
	
	//DOF contains local-global index 
	// no shape functions now.
	
//	public void assignTo(Element e) {
//		e.clearAllDOF();
//		if(nTotalNodes == -1 || nDOF_p == -1) {
//			FutureyeException ex = new FutureyeException("Call initDOFIndex() first!");
//			ex.printStackTrace();
//			System.exit(-1);
//		}
//		//单元结点数
//		int nNode = e.nodes.size();
//		//Assign shape function to DOF
//		for(int j=1;j<=nNode;j++) {
//			//Asign shape function to DOF
//			DOF dof_u1 = new DOF(
//					j,//Local DOF index
//					//Global DOF index, take global node index
//					e.nodes.at(j).globalIndex,
//					null
//					         );
//			dof_u1.setVVFComponent(1);
//			DOF dof_u2 = new DOF(
//					nNode+j,//Local DOF index
//					//Global DOF index, take this.nTotalNodes + global node index
//					this.nTotalNodes+e.nodes.at(j).globalIndex,
//					null
//					         );
//			dof_u2.setVVFComponent(2);
//			e.addNodeDOF(j, dof_u1);
//			//e.addNodeDOF(j, dof_u2); //???bug???
//			e.addNodeDOF(nNode+j, dof_u2);
//		}
//		
//		//Assign shape function to DOF
//		DOF dof = new DOF(
//					2*nNode+1, //Local DOF index
//					//this.nTotalNodes*2+nDOF_p, //Global DOF index for Pressure
//					this.nTotalNodes*2+this.nDOF_p, //Global DOF index for Pressure
//					shapeFun[2*nNode] //Shape function 
//					);
//		this.nDOF_p++;
//		dof.setVVFComponent(3);	
//		e.addVolumeDOF(dof);
//	}

	@Override
	public VecFiniteElement getBoundaryFE() {
		return new FELinearV_ConstantPLine2D();
	}

	@Override
	public boolean isDOFCoupled(int idx1, int idx2) {
		if(idx1 == 8 || idx2 == 8)
			return true;
		if(idx1 <= 3 && idx2 >= 4)
			return false;
		if(idx2 <= 3 && idx1 >= 4)
			return false;
		return true;
	}
}
