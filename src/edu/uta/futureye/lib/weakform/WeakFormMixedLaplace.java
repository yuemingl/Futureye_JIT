package edu.uta.futureye.lib.weakform;

import edu.uta.futureye.algebra.intf.BlockMatrix;
import edu.uta.futureye.algebra.intf.BlockVector;
import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.core.DOF;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.function.intf.VectorShapeFunction;
import edu.uta.futureye.function.operator.FMath;
import edu.uta.futureye.function.operator.FOIntegrate;
import edu.uta.futureye.util.Utils;
import edu.uta.futureye.util.container.DOFList;

/**
 * Problem:
 * 
 *   div \mathbf{p} + f = 0
 *   \mathbf{p} = \nabla u
 * 
 * Weak Form: 
 *  seek \mathbf{p} \in H_{g,N}(div,\Omega) and u \in L^{2}(\Omega)
 *  such that
 *  (p,q)_{\Omega} + (u,\div{q})_{\Omega} = u_D*(q,\mu)_{\Gamma_{D}},  
 *  							for all q \in H_{0,N}(div,\Omega)
 *  (v,\div{p})_{\Omega} = -(v,f)_{\Omega}, 
 *  							for all v \in L^{2}(\Omega)
 * 
 * Algebra System:
 * 
 * ( B  C )(p)   (b0)
 * ( C' 0 )(u) = (bf)
 * 
 * @author liuyueming
 *
 */
public class WeakFormMixedLaplace extends AbstractVectorWeakform {
	protected Function g_f = null;
	protected Function g_k = null;
	//protected Function g_c = null;
	//protected Function g_g = null;
	//protected Function g_d = null;
	
	public void setF(Function f) {
		this.g_f = f;
	}
	
	//Robin:  d*u + k*u_n = g
	public void setParam(Function k,Function c,Function g,Function d) {
		this.g_k = k;
		//this.g_c = c;
		//this.g_g = g;
		//this.g_d = d;
	}
	
	@Override
	public void assembleElement(Element e, 
			Matrix globalStiff, Vector globalLoad){
		
		DOFList edgeDOFs = e.getAllEdgeDOFList();
		//获取单元体对应的自由度列表
		DOFList eleDOFs = e.getVolumeDOFList();
		int nEdgeDOF = edgeDOFs.size();
		int nElementDOF = eleDOFs.size();

		BlockMatrix blockMat = (BlockMatrix)globalStiff;
		BlockVector blockVec = (BlockVector)globalLoad;
		
//不需要使用单独的每个分块子矩阵，直接使用整个块矩阵就可以了，
//合成过程由于整个块矩阵的特性，会自动将相应的元素放入对应的子块中。		
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
				
				//B = (p,q)_{\Omega}
				Function integrandB = null;
				integrandB = vecU.dot(vecV);
				//单元上数值积分
				if(e.vertices().size() == 3) {
					double val = FOIntegrate.intOnTriangleRefElement(
							integrandB.M(e.getJacobin()),4);
					blockMat.add(dofU.getGlobalIndex(), dofV.getGlobalIndex(), val);
				}
			}
			//面自由度循环（2D单元）
			for(int k=1;k<=nElementDOF;k++) {
				DOF dofE = eleDOFs.at(k);
				//C = (u,\div{q})_{\Omega}
				Function integrandC = null;
				integrandC = FMath.div(vecV);
				double val = FOIntegrate.intOnTriangleRefElement(
						integrandC.M(e.getJacobin()),4);
				blockMat.add(dofV.getGlobalIndex(), dofE.getGlobalIndex(), val);
				//C' = (v,\div{p})_{\Omega}
				blockMat.add(dofE.getGlobalIndex(), dofV.getGlobalIndex(), val);
			}
			
			//b0 = 0 //Dirichlet条件如何处理？
		}
		
		for(int k=1;k<=nElementDOF;k++) {
			DOF dofE = eleDOFs.at(k);
			//ShapeFunction sf = dofE.getSF(); //分片常数元，在积分项中系数是1
			Function integrand = Utils.interpolateFunctionOnElement(g_f, e);
			//bf = -(v,f)_{\Omega}
			integrand = FC.c(-1.0).M(integrand);
			double val = 0.0;
			if(e.vertices().size() == 3) {
				val = FOIntegrate.intOnTriangleRefElement(
						integrand.M(e.getJacobin()),4
					);
			} else if (e.vertices().size() == 4) {
				val = FOIntegrate.intOnRectangleRefElement(
						integrand.M(e.getJacobin()),2 //TODO
						);
			}
			blockVec.add(dofE.getGlobalIndex(), val);
		}
	}
}
