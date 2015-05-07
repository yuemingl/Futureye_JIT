package edu.uta.futureye.lib.weakform;

import java.util.HashMap;
import java.util.Map;

import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.core.DOF;
import edu.uta.futureye.core.DOFOrder;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.NodeType;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.function.intf.ScalarShapeFunction;
import edu.uta.futureye.function.operator.FMath;
import edu.uta.futureye.function.operator.FOIntegrate;
import edu.uta.futureye.util.Utils;
import edu.uta.futureye.util.container.DOFList;
import edu.uta.futureye.util.container.ElementList;

/**
 * <blockquote><pre>
 * Solve
 *   -k*Laplace(u) + c*u = f, in \Omega
 *   u = u0,                  on \Gamma1
 *   d*u + k*u_n = g,         on \Gamma2
 *=>
 *   A(u, v) = (f, v)
 * 
 * where
 *   A(u, v) = (k*u_x, v_x) + (k*u_y, v_y) + (k*u_z, v_z) - (g-d*u,v)_\Gamma2 + (c*u, v)
 *=>
 *   A(u, v) = (k*u_x, v_x) + (k*u_y, v_y) + (k*u_z, v_z) + (d*u-g,v)_\Gamma2 + (c*u, v)
 *
 *   \Gamma1: Dirichlet boundary of \Omega
 *   \Gamma2: Neumann(Robin) boundary of \Omega
 *   u_n: \frac{\pratial{u}}{\partial{n}}
 *   n: unit norm vector of \Omega
 *   k = k(x,y)
 *   c = c(x,y)
 *   d = d(x,y)
 *   g = g(x,y)
 * </blockquote></pre>  
 * 
 * @author liuyueming
 *
 */
public class WeakFormLaplace3D extends AbstractScalarWeakForm {
	protected MathFunc g_f = null;
	protected MathFunc g_k = null;
	protected MathFunc g_c = null;
	protected MathFunc g_g = null;
	protected MathFunc g_d = null;

	public void setF(MathFunc f) {
		this.g_f = f;
	}
	
	//Robin:  d*u + k*u_n= g (自然边界：d==k, g=0)
	public void setParam(MathFunc k,MathFunc c,MathFunc g,MathFunc d) {
		this.g_k = k;
		this.g_c = c;
		this.g_g = g;
		this.g_d = d;
	}

	@Override
	public MathFunc leftHandSide(Element e, ItemType itemType) {
		if(itemType==ItemType.Domain)  {
			MathFunc integrand = null;
			if(g_k == null) {
				integrand = FMath.sum(
						u.diff("x").M(v.diff("x")),
						u.diff("y").M(v.diff("y")),
						u.diff("z").M(v.diff("z"))
					);
			} else {
				MathFunc fk = Utils.interpolateOnElement(g_k,e);
				MathFunc fc = Utils.interpolateOnElement(g_c,e);
				integrand = fk.M(FMath.sum(
								u.diff("x").M(v.diff("x")),
								u.diff("y").M(v.diff("y")),
								u.diff("z").M(v.diff("z"))
							)).A(fc.M(u.M(v)));
			}
			return integrand;
		} else if(itemType==ItemType.Border) {//Neumann border integration on LHS
			if(g_d != null) {
				Element be = e;
				MathFunc fd = Utils.interpolateOnElement(g_d, be);
				MathFunc borderIntegrand = fd.M(u.M(v));
				return borderIntegrand;
			}
		}
		return null;
	}

	@Override
	public MathFunc rightHandSide(Element e, ItemType itemType) {
		if(itemType==ItemType.Domain)  {
			MathFunc ff = Utils.interpolateOnElement(g_f, e);
			MathFunc integrand = ff.M(v);
			return integrand;
		} else if(itemType==ItemType.Border) {//Neumann border type
			if(g_g != null) {
				Element be = e;
				MathFunc fq = Utils.interpolateOnElement(g_g, be);
				MathFunc borderIntegrand = fq.M(v);
				return borderIntegrand;
			}
		}
		return  null;
	}
	
	
	/**
	 * Optimized for fast assemble, 10% speedup
	 */
	@Override
	public void assembleElement(Element e, 
			Matrix globalStiff, Vector globalLoad) {
		
		DOFList DOFs = e.getAllDOFList(DOFOrder.NEFV);
		int nDOFs = DOFs.size();
		
		//Update Jacobin on e
		e.updateJacobinLinear3D();
		
		//形函数计算需要和单元关联，并提前计算导数
		Map<Integer, MathFunc> mapShape_x = new HashMap<Integer, MathFunc>();
		Map<Integer, MathFunc> mapShape_y = new HashMap<Integer, MathFunc>();
		Map<Integer, MathFunc> mapShape_z = new HashMap<Integer, MathFunc>();
		for(int i=1;i<=nDOFs;i++) {
			DOF dof = DOFs.at(i);
			ScalarShapeFunction sf = dof.getSSF();
			dof.getSSF().assignElement(e);
			mapShape_x.put(dof.getLocalIndex(), sf.diff("x"));
			mapShape_y.put(dof.getLocalIndex(), sf.diff("y"));
			mapShape_z.put(dof.getLocalIndex(), sf.diff("z"));
		}

		MathFunc fk = null;
		if(g_k != null) fk = Utils.interpolateOnElement(g_k,e);
		MathFunc fc = null;
		if(g_c != null) fc = Utils.interpolateOnElement(g_c,e);

		//所有自由度双循环
		for(int i=1;i<=nDOFs;i++) {
			DOF dofI = DOFs.at(i);
			ScalarShapeFunction sfI = dofI.getSSF();
			int nLocalRow = dofI.getLocalIndex();
			int nGlobalRow = dofI.getGlobalIndex();
			for(int j=1;j<=nDOFs;j++) {
				DOF dofJ = DOFs.at(j);
				int nLocalCol = dofJ.getLocalIndex();
				int nGlobalCol = dofJ.getGlobalIndex();
				
				//Integrand part of Weak Form on element e
				MathFunc integrand = null;
				if(g_k == null) {
					integrand = FMath.sum(
						mapShape_x.get(nLocalRow).M(mapShape_x.get(nLocalCol)),
						mapShape_y.get(nLocalRow).M(mapShape_y.get(nLocalCol)),
						mapShape_z.get(nLocalRow).M(mapShape_z.get(nLocalCol))
						);
				} else {
					integrand = fk.M(FMath.sum(
									mapShape_x.get(nLocalRow).M(mapShape_x.get(nLocalCol)),
									mapShape_y.get(nLocalRow).M(mapShape_y.get(nLocalCol)),
									mapShape_z.get(nLocalRow).M(mapShape_z.get(nLocalCol))
									)
								).A(
									fc.M(dofI.getSSF().M(dofJ.getSSF()))
								);
				}
				//Numerical integration on element e
				double lhsVal = 0.0;
				if (e.vertices().size() == 4) {
					lhsVal = FOIntegrate.intOnTetrahedraRefElement(
							integrand.M(e.getJacobin()),2
							);
				}
				globalStiff.add(nGlobalRow, nGlobalCol, lhsVal);
			}
			//Load vector
			MathFunc ff = Utils.interpolateOnElement(g_f, e);
			MathFunc integrand = ff.M(sfI);
			double rhsVal = 0.0;
			if (e.vertices().size() == 4) {
				rhsVal = FOIntegrate.intOnTetrahedraRefElement(
						integrand.M(e.getJacobin()),2
						);
			}
			globalLoad.add(nGlobalRow, rhsVal);
		}
		
		//Robin:  d*u + k*u_n= q (自然边界：d==k, q=0)
		if(g_d != null && e.isBorderElement()) {
			ElementList beList = e.getBorderElements();
			for(int n=1;n<=beList.size();n++) {
				Element be = beList.at(n);
				
				MathFunc fd = null;
				if(g_d != null) fd = Utils.interpolateOnElement(g_d, be);
				
				//Check node type
				NodeType nodeType = be.getBorderNodeType();
				if(nodeType == NodeType.Neumann || nodeType == NodeType.Robin) {
					
					//Update Jacobin on be
					be.updateJacobinLinear2D();
					DOFList beDOFs = be.getAllDOFList(DOFOrder.NEFV);
					int nBeDOF = beDOFs.size();
					
					//形函数计算需要和单元关联
					for(int i=1;i<=nBeDOF;i++) {
						beDOFs.at(i).getSSF().assignElement(be);
					}
					
					//所有自由度双循环
					for(int i=1;i<=nBeDOF;i++) {
						DOF dofI = beDOFs.at(i);
						ScalarShapeFunction sfI = dofI.getSSF();
						int nGlobalRow = dofI.getGlobalIndex();
						for(int j=1;j<=nBeDOF;j++) {
							DOF dofJ = beDOFs.at(j);
							ScalarShapeFunction sfJ = dofJ.getSSF();
							int nGlobalCol = dofJ.getGlobalIndex();
							//Stiff matrix for border
							MathFunc borderIntegrand = fd.M(sfI.M(sfJ));
							//Numerical integrate the border 'be' of element 'e'
							double lhsBrVal = 0.0;
							if(be.vertices().size() == 3) {
								lhsBrVal = FOIntegrate.intOnTriangleRefElement(
											borderIntegrand.M(be.getJacobin()),5
										);
							} else if (be.vertices().size() == 4) {
								lhsBrVal = FOIntegrate.intOnRectangleRefElement(
											borderIntegrand.M(be.getJacobin()),2 //TODO
										);
							}
							globalStiff.add(nGlobalRow, nGlobalCol, lhsBrVal);
						}
						//Load vector for border
						if(g_g != null) {
							MathFunc fq = Utils.interpolateOnElement(g_g, be);
							MathFunc borderIntegrand = fq.M(v);
							double rhsBrVal = 0.0;
							if(be.vertices().size() == 3) {
								rhsBrVal = FOIntegrate.intOnTriangleRefElement(
											borderIntegrand.M(be.getJacobin()),5
										);
							} else if (be.vertices().size() == 4) {
								rhsBrVal = FOIntegrate.intOnRectangleRefElement(
											borderIntegrand.M(be.getJacobin()),2 //TODO
										);
							}
							globalLoad.add(nGlobalRow, rhsBrVal);	
						}
					}
				}
			}
		}
	}
}