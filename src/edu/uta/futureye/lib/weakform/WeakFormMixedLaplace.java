package edu.uta.futureye.lib.weakform;

import edu.uta.futureye.algebra.intf.BlockMatrix;
import edu.uta.futureye.algebra.intf.BlockVector;
import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.core.DOF;
import edu.uta.futureye.core.DOFOrder;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.intf.WeakForm;
import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.function.intf.ShapeFunction;
import edu.uta.futureye.function.intf.VectorShapeFunction;
import edu.uta.futureye.function.operator.FOBasic;
import edu.uta.futureye.function.operator.FOIntegrate;
import edu.uta.futureye.function.operator.FOVector;
import edu.uta.futureye.util.Utils;
import edu.uta.futureye.util.list.DOFList;

/**
 * Problem:
 * 
 *   div \mathbf{p} + f = 0
 *   \mathbf{p} = \Nabla u
 * 
 * Weak Form: 
 *  seek q \in H_{g,N}(div,\Omega) and u \in L^{2}(\Omega)
 *  such that
 *  (p,q)_{\Omega} + (u,\div{q})_{\Omega} = u_D*(q,\mu)_{\Gamma_{D}},  for all q \in H_{0,N}(div,\Omega)
 *  (v,\div{p})_{\Omega} = -(v,f)_{\Omega}, for all v \in L^{2}(\Omega)
 * 
 * Algebra System:
 * 
 * ( B  C )   (b0)
 * ( C' 0 ) = (bf)
 * 
 * @author liuyueming
 *
 */
public class WeakFormMixedLaplace implements WeakForm {
	protected ShapeFunction u = null;
	protected ShapeFunction v = null;
	
	protected int uDOFLocalIndex;
	protected int vDOFLocalIndex;
	

	protected Function g_f = null;
	protected Function g_k = null;
	protected Function g_c = null;
	protected Function g_q = null;
	protected Function g_d = null;
	
	@Override
	public void assembleElement(Element e, 
			Matrix globalStiff, Vector globalLoad){
		
		DOFList edgeDOFs = e.getAllEdgeDOFList(DOFOrder.NEFV);
		DOFList eleDOFs = e.getElementDOFList();
		int nEdgeDOF = edgeDOFs.size();
		int nElementDOF = eleDOFs.size();

		BlockMatrix blockMat = (BlockMatrix)globalStiff;
		BlockVector blockVec = (BlockVector)globalLoad;
//不需要使用分块矩阵，直接使用整个矩阵就可以了，合成过程由于块矩阵的性质
//会自动将相应的元素放入对应的块中。		
//		Matrix m11 = blockMat.getBlock(1, 1);
//		Matrix m12 = blockMat.getBlock(1, 2);
//		Matrix m21 = blockMat.getBlock(2, 1);
//		//Matrix m22 = blockMat.getBlock(2, 2);
		
		e.updateJacobinLinear2D();
		for(int i=1;i<=nEdgeDOF;i++) {
			edgeDOFs.at(i).getVSF().asignElement(e);
		}
		
		//边自由度双循环
		for(int j=1;j<=nEdgeDOF;j++) {
			DOF dofV = edgeDOFs.at(j);
			VectorShapeFunction vecV = dofV.getVSF();
			for(int i=1;i<=nEdgeDOF;i++) {
				DOF dofU = edgeDOFs.at(i);
				VectorShapeFunction vecU = dofU.getVSF();
				
				//B: (p,q)_{\Omega}
				Function integrandB = null;
				integrandB = vecU.dot(vecV);
				//单元上数值积分
				Function integralB = null;
				if(e.vertices().size() == 3) {
					integralB = FOIntegrate.intOnTriangleRefElement(
							FOBasic.Mult(integrandB, e.getJacobin()),5
					);
					double val = integralB.value(null);
					blockMat.add(dofU.getGlobalIndex(), dofV.getGlobalIndex(), val);
				}
			}
			//面自由度循环（2D单元）
			for(int k=1;k<=nElementDOF;k++) {
				DOF dofE = eleDOFs.at(k);
				//C: (u,\div{q})_{\Omega}
				Function integrandC = null;
				integrandC = FOVector.Div(vecV);
				Function integralC = FOIntegrate.intOnTriangleRefElement(
						FOBasic.Mult(integrandC, e.getJacobin()),5
						);
				double val = integralC.value(null);
				blockMat.add(dofV.getGlobalIndex(), dofE.getGlobalIndex(), val);
				//C': (v,\div{p})_{\Omega}
				blockMat.add(dofE.getGlobalIndex(), dofV.getGlobalIndex(), val);
			}
			
			//b0 Dirichlet条件？
		}
		
		for(int k=1;k<=nElementDOF;k++) {
			DOF dofE = eleDOFs.at(k);
			//ShapeFunction sf = dofE.getSF(); //分片常数元，在积分项中系数是1
			Function integrand = Utils.interplateFunctionOnElement(g_f, e);
			//bf: -(v,f)_{\Omega}
			integrand = FOBasic.Mult(new FC(-1.0), integrand);
			Function integral = null;
			if(e.vertices().size() == 3) {
				integral = FOIntegrate.intOnTriangleRefElement(
						FOBasic.Mult(integrand, e.getJacobin()),5
					);
			} else if (e.vertices().size() == 4) {
				integral = FOIntegrate.intOnRectangleRefElement(
						FOBasic.Mult(integrand, e.getJacobin()),2 //TODO
						);
			}
			double val = integral.value(null);
			blockVec.add(dofE.getGlobalIndex(), val);
		}
	}
	
	@Override
	public Function leftHandSide(Element e, ItemType itemType) {
		return null;
	}

	@Override
	public Function rightHandSide(Element e, ItemType itemType) {
		return null;
	}

	@Override
	public void setShapeFunction(ShapeFunction trial, int trialDofLocalIndex,
			ShapeFunction test, int testDofLocalIndex) {
	}
	
	public void setF(Function f) {
		this.g_f = f;
	}
	
	//Robin:  d*u + k*u_n= q (自然边界：d==k, q=0)
	public void setParam(Function k,Function c,Function q,Function d) {
		this.g_k = k;
		this.g_c = c;
		this.g_q = q;
		this.g_d = d;
	}

	@Override
	public double integrate(Element e, Function fun) {
		return 0;
	}
}
