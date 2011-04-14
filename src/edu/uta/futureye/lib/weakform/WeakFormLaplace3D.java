package edu.uta.futureye.lib.weakform;

import java.util.HashMap;
import java.util.Map;

import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.core.DOF;
import edu.uta.futureye.core.DOFOrder;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.NodeType;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.function.intf.ScalarShapeFunction;
import edu.uta.futureye.function.operator.FMath;
import edu.uta.futureye.function.operator.FOIntegrate;
import edu.uta.futureye.util.Utils;
import edu.uta.futureye.util.container.DOFList;
import edu.uta.futureye.util.container.ElementList;

/**
 * Solve
 *   -k*Laplace(u) + c*u = f, in \Omega
 *   u = u0,                  on \Gamma1
 *   d*u + k*u_n = q,         on \Gamma2
 *=>
 *   A(u, v) = (f, v)
 * 
 * where
 *   A(u, v) = (k*u_x, v_x) + (k*u_y, v_y) + (k*u_z, v_z) - (q-d*u,v)_\Gamma2 + (c*u, v)
 *=>
 *   A(u, v) = (k*u_x, v_x) + (k*u_y, v_y) + (k*u_z, v_z) + (d*u-q,v)_\Gamma2 + (c*u, v)
 *
 *   \Gamma1: Dirichlet boundary of \Omega
 *   \Gamma2: Neumann(Robin) boundary of \Omega
 *   u_n: \frac{\pratial{u}}{\partial{n}}
 *   n: unit norm vector of \Omega
 *   k = k(x,y)
 *   c = c(x,y)
 *   d = d(x,y)
 *   q = q(x,y)
 *   
 * @author liuyueming
 *
 */
public class WeakFormLaplace3D extends AbstractScalarWeakForm {
	protected Function g_f = null;
	protected Function g_k = null;
	protected Function g_c = null;
	protected Function g_q = null;
	protected Function g_d = null;

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
	public Function leftHandSide(Element e, ItemType itemType) {
		if(itemType==ItemType.Domain)  {
			Function integrand = null;
			if(g_k == null) {
				integrand = FMath.sum(
						u._d("x").M(v._d("x")),
						u._d("y").M(v._d("y")),
						u._d("z").M(v._d("z"))
					);
			} else {
				Function fk = Utils.interplateFunctionOnElement(g_k,e);
				Function fc = Utils.interplateFunctionOnElement(g_c,e);
				integrand = fk.M(FMath.sum(
								u._d("x").M(v._d("x")),
								u._d("y").M(v._d("y")),
								u._d("z").M(v._d("z"))
							)).A(fc.M(u.M(v)));
			}
			return integrand;
		} else if(itemType==ItemType.Border) {//Neumann border integration on LHS
			if(g_d != null) {
				Element be = e;
				Function fd = Utils.interplateFunctionOnElement(g_d, be);
				Function borderIntegrand = fd.M(u.M(v));
				return borderIntegrand;
			}
		}
		return null;
	}

	@Override
	public Function rightHandSide(Element e, ItemType itemType) {
		if(itemType==ItemType.Domain)  {
			Function ff = Utils.interplateFunctionOnElement(g_f, e);
			Function integrand = ff.M(v);
			return integrand;
		} else if(itemType==ItemType.Border) {//Neumann border type
			Element be = e;
			Function fq = Utils.interplateFunctionOnElement(g_q, be);
			Function borderIntegrand = fq.M(v);
			return borderIntegrand;
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
		Map<Integer, Function> mapShape_x = new HashMap<Integer, Function>();
		Map<Integer, Function> mapShape_y = new HashMap<Integer, Function>();
		Map<Integer, Function> mapShape_z = new HashMap<Integer, Function>();
		for(int i=1;i<=nDOFs;i++) {
			DOF dof = DOFs.at(i);
			ScalarShapeFunction sf = dof.getSSF();
			dof.getSSF().asignElement(e);
			mapShape_x.put(dof.getLocalIndex(), sf._d("x"));
			mapShape_y.put(dof.getLocalIndex(), sf._d("y"));
			mapShape_z.put(dof.getLocalIndex(), sf._d("z"));
		}

		Function fk = null;
		if(g_k != null) fk = Utils.interplateFunctionOnElement(g_k,e);
		Function fc = null;
		if(g_c != null) fc = Utils.interplateFunctionOnElement(g_c,e);

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
				Function integrand = null;
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
				Function integral = null;
				if (e.vertices().size() == 4) {
					integral = FOIntegrate.intOnTetrahedraRefElement(
							integrand.M(e.getJacobin()),2
							);
				}
				double lhsVal = integral.value(null);
				globalStiff.add(nGlobalRow, nGlobalCol, lhsVal);
			}
			//Load vector
			Function ff = Utils.interplateFunctionOnElement(g_f, e);
			Function integrand = ff.M(sfI);
			Function integral = null;
			if (e.vertices().size() == 4) {
				integral = FOIntegrate.intOnTetrahedraRefElement(
						integrand.M(e.getJacobin()),2
						);
			}
			double rhsVal = integral.value(null);
			globalLoad.add(nGlobalRow, rhsVal);
		}
		
		//Robin:  d*u + k*u_n= q (自然边界：d==k, q=0)
		if(g_d != null && e.isBorderElement()) {
			ElementList beList = e.getBorderElements();
			for(int n=1;n<=beList.size();n++) {
				Element be = beList.at(n);
				
				Function fd = null;
				if(g_d != null) fd = Utils.interplateFunctionOnElement(g_d, be);
				
				//Check node type
				NodeType nodeType = be.getBorderNodeType();
				if(nodeType == NodeType.Neumann || nodeType == NodeType.Robin) {
					
					//Update Jacobin on be
					be.updateJacobinLinear2D();
					DOFList beDOFs = be.getAllDOFList(DOFOrder.NEFV);
					int nBeDOF = beDOFs.size();
					
					//形函数计算需要和单元关联
					for(int i=1;i<=nBeDOF;i++) {
						beDOFs.at(i).getSSF().asignElement(be);
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
							Function borderIntegrand = fd.M(sfI.M(sfJ));
							//Numerical integrate the border 'be' of element 'e'
							Function borderIntgral = null;
							if(be.vertices().size() == 3) {
								borderIntgral = FOIntegrate.intOnTriangleRefElement(
											borderIntegrand.M(be.getJacobin()),5
										);
							} else if (be.vertices().size() == 4) {
								borderIntgral = FOIntegrate.intOnRectangleRefElement(
											borderIntegrand.M(be.getJacobin()),2 //TODO
										);
							}
							double lhsBrVal = borderIntgral.value(null);
							globalStiff.add(nGlobalRow, nGlobalCol, lhsBrVal);
						}
						//Load vector for border
						if(g_q != null) {
							Function fq = Utils.interplateFunctionOnElement(g_q, be);
							Function borderIntegrand = fq.M(v);
							Function borderIntgral = null;
							if(be.vertices().size() == 3) {
								borderIntgral = FOIntegrate.intOnTriangleRefElement(
											borderIntegrand.M(be.getJacobin()),5
										);
							} else if (be.vertices().size() == 4) {
								borderIntgral = FOIntegrate.intOnRectangleRefElement(
											borderIntegrand.M(be.getJacobin()),2 //TODO
										);
							}
							double rhsBrVal = borderIntgral.value(null);
							globalLoad.add(nGlobalRow, rhsBrVal);	
						}
					}
				}
			}
		}
	}
}