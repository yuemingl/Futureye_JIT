package edu.uta.futureye.lib.weakform;

import static symjava.math.SymMath.dot;
import static symjava.math.SymMath.grad;
import static symjava.symbolic.Symbol.r;
import static symjava.symbolic.Symbol.s;
import static symjava.symbolic.Symbol.u;
import static symjava.symbolic.Symbol.v;
import static symjava.symbolic.Symbol.x;
import static symjava.symbolic.Symbol.y;

import java.util.HashMap;
import java.util.Map;

import symjava.examples.fem.UnitRightTriangle;
import symjava.math.Transformation;
import symjava.matrix.SymMatrix;
import symjava.numeric.NumInt;
import symjava.relational.Eq;
import symjava.symbolic.Expr;
import symjava.symbolic.Func;
import symjava.symbolic.Int;
import symjava.symbolic.SymConst;
import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.core.DOF;
import edu.uta.futureye.core.DOFOrder;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.NodeType;
import edu.uta.futureye.function.FEMFunc;
import edu.uta.futureye.function.intf.ScalarShapeFunction;
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
 *   A(u, v) = (k*u_x, v_x) + (k*u_y, v_y) - (g-d*u,v)_\Gamma2 + (c*u, v)
 *=>
 *   A(u, v) = (k*u_x, v_x) + (k*u_y, v_y) + (d*u-g,v)_\Gamma2 + (c*u, v)
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
public class SymWeakFormLaplace2D extends AbstractScalarWeakForm {
//	protected Expr g_f = null;
//	protected Expr g_k = null;
//	protected Expr g_c = null;
//	protected Expr g_g = null; 
//	protected Expr g_d = null;
	
	protected Expr ff = null;
	protected Expr fk = null;
	protected Expr fc = null;
	protected Expr fg = null;
	protected Expr fd = null;
	
	protected FEMFunc lhs = null;
	protected FEMFunc rhs = null;
	
//	NumInt lhsNInt[][] = new NumInt[shapeFuns.length][shapeFuns.length];
//	NumInt rhsNInt[] = new NumInt[shapeFuns.length];
//	NumInt lhsNInt[][] = new NumInt[dof.length][dof.length];
//	NumInt rhsNInt[] = new NumInt[dof.length];	
	FEMFunc lhsNInt[][] = new FEMFunc[3][3];
	FEMFunc rhsNInt[] = new FEMFunc[3];	
	
	public void setF(Expr f) {
		this.ff = f;
	}
	
	//Robin:  d*u + k*u_n = g 
	//2011-08-02
	//E.g.1 Nature boundary condition: u_n + u = 0  =>  d=k, q=0
	//E.g.2 u_n = g                                 =>  d=0, g=k*g  
	public void setParam(Expr k,Expr c,Expr g,Expr d) {
		this.fk = k;
		this.fc = c;
		this.fg = g;
		this.fd = d;
	}
	
	@Override 
	public void preProcess(Element e) {
//		if(e.dim() == 2) {
//			if(g_k != null) fk = Utils.interpolateOnElement(g_k, e);
//			if(g_c != null) fc = Utils.interpolateOnElement(g_c, e);
//			if(g_f != null) ff = Utils.interpolateOnElement(g_f, e);
//		} else if(e.dim() == 1) {
//			if(g_d != null) fd = Utils.interpolateOnElement(g_d, e);
//			if(g_g != null) fg = Utils.interpolateOnElement(g_g, e);
//		}
	}

	public void generateOrLoadBytecode() {
		Func u = new Func("u", x, y);
		Func v = new Func("v", x, y);
		Expr integrandLHS = fk*dot(grad(u), grad(v))+fc*u*v;
		Expr integrandRHS = ff*v;
		
		//No need: e.getTransFormation
		//e.getJacobian
		//Create coordinate transformation on template element
		SymConst x1 = new SymConst("x1");
		SymConst x2 = new SymConst("x2");
		SymConst x3 = new SymConst("x3");
		SymConst y1 = new SymConst("y1");
		SymConst y2 = new SymConst("y2");
		SymConst y3 = new SymConst("y3");
		Transformation trans = new Transformation(
				new Eq(x, x1*r+x2*s+x3*(1-r-s)),
				new Eq(y, y1*r+y2*s+y3*(1-r-s))
				);
		// jac = (xr xs)
		//       (yr ys)
		SymMatrix jacMat = trans.getJacobianMatrix();
		System.out.println(jacMat);
		System.out.println();
		
		
		//e.getShapeFunctions
		//e.getDOFAllList
		//Shape functions
		Func N1 = new Func("R", x, y);
		Func N2 = new Func("S", x, y);
		Func N3 = new Func("T", 1 - N1 - N2);
		Func[] shapeFuns = {N1, N2, N3};
		
		Expr jac = trans.getJacobian();
		Expr rx =  jacMat[1][1]/jac; //rx =  ys/jac
		Expr ry = -jacMat[0][1]/jac; //ry = -xs/jac //bugfix missed a minus sign!
		Expr sx = -jacMat[1][0]/jac; //sx = -yr/jac
		Expr sy =  jacMat[0][0]/jac; //sy =  xr/jac
		System.out.println(jac);
		System.out.println(rx);
		System.out.println(ry);
		System.out.println(sx);
		System.out.println(sy);

		//This is the reference element: a unit right triangle
		UnitRightTriangle tri = new UnitRightTriangle("Tri", r, s);
		
		Int lhsInt[][] = new Int[shapeFuns.length][shapeFuns.length];
		Int rhsInt[] = new Int[shapeFuns.length];
		for(int i=0; i<shapeFuns.length; i++) {
			Func V = shapeFuns[i]; //test
			for(int j=0; j<shapeFuns.length; j++) {
				Func U = shapeFuns[j]; //trial
				
				//Weak form for the left hand side of the PDE
				Expr lhs = integrandLHS.subs(u, U).subs(v, V);
				System.out.println(lhs);
				System.out.println();
				
				//Replace the derivatives with it's concrete expression
				lhs = lhs
					.subs(N1.diff(x), rx).subs(N1.diff(y), ry)
					.subs(N2.diff(x), sx).subs(N2.diff(y), sy)
					.subs(N1, r).subs(N2, s)
					.subs(x, trans.eqs[0].rhs)
					.subs(y, trans.eqs[1].rhs);
				System.out.println(lhs);
				System.out.println();
				
				//Define the integration on the reference domain
				//lhsInt[i][j] = new Int(new Func("",lhs*jac,new Expr[]{r,s}), tri);
				lhsInt[i][j] = new Int(new Func(
						String.format("lhs%d%d",i,j), lhs*jac,new Expr[]{r,s}), tri);
				System.out.println(lhsInt[i][j]);
				System.out.println();
			}
			
			
			//Weak form for the right hand side of the PDE
			Expr rhs = integrandRHS.subs(v, V)
					.subs(N1, r).subs(N2, s)
					.subs(x, trans.eqs[0].rhs)
					.subs(y, trans.eqs[1].rhs);
			//System.out.println(rhs);
			System.out.println();
			rhsInt[i] = new Int(new Func(
					String.format("rhs%d",i),rhs*jac,new Expr[]{r,s}), tri);
		}
		
		//Generate bytecode for the integration
		//You can save the class to some place for later use
		for(int i=0; i<shapeFuns.length; i++) {
			for(int j=0; j<shapeFuns.length; j++) {
				lhsNInt[i][j] = new FEMFunc(lhsInt[i][j]);
			}
			rhsNInt[i] = new FEMFunc(rhsInt[i]);
		}
		
	}
	
	@Override
	public FEMFunc leftHandSide(Element e, ItemType itemType) {
		if(itemType==ItemType.Domain)  {
			//Integrand part of Weak Form on element e
			int idxI = this.trialDOF.getLocalIndex();
			int idxJ = this.testDOF.getLocalIndex();
			return this.lhsNInt[idxJ-1][idxI-1];
//Old code			
//			if(g_k == null) {
//				integrand = u._d("x").M(v._d("x")) .A (u._d("y").M(v._d("y")));
//			} else {
//				integrand = fk.M(
//								u._d("x").M(v._d("x")) .A (u._d("y").M(v._d("y")))
//							).A(
//								fc.M(u.M(v))
//							);
//			}
//			return integrand;
		}
//		else if(itemType==ItemType.Border) {//Neumann border integration on LHS
//			if(g_d != null) {
//				FEMFunc borderIntegrand = fd.M(u.M(v));
//				return borderIntegrand;
//			}
//		}
		return null;
	}

	@Override
	public FEMFunc rightHandSide(Element e, ItemType itemType) {
		if(itemType==ItemType.Domain)  {
			int idxJ = this.testDOF.getLocalIndex();
			return rhsNInt[idxJ-1];
			
//			FEMFunc integrand = ff.M(v);
//			return integrand;
		} 
//		else if(itemType==ItemType.Border) {
//			if(g_g != null) {
//				FEMFunc borderIntegrand = fg.M(v);
//				return borderIntegrand;
//			}
//		}
		return null;	
	}

//	/**
//	 * Optimized for fast assemble, 10% speedup
//	 */
//	@Override
//	public void assembleElement(Element e, 
//		Matrix globalStiff,	Vector globalLoad) {
//	
//		DOFList DOFs = e.getAllDOFList(DOFOrder.NEFV);
//		int nDOFs = DOFs.size();
//		
////		//Update Jacobin on e
////		e.updateJacobinLinear2D();
//		
////		//形函数计算需要和单元关联，并提前计算导数
////		Map<Integer, FEMFunc> mapShape_x = new HashMap<Integer, FEMFunc>();
////		Map<Integer, FEMFunc> mapShape_y = new HashMap<Integer, FEMFunc>();
////		for(int i=1;i<=nDOFs;i++) {
////			DOF dof = DOFs.at(i);
////			ScalarShapeFunction sf = dof.getSSF();
////			dof.getSSF().assignElement(e);
////			mapShape_x.put(dof.getLocalIndex(), sf._d("x"));
////			mapShape_y.put(dof.getLocalIndex(), sf._d("y"));
////		}
//
////		FEMFunc fk = null;
////		if(g_k != null) fk = Utils.interpolateOnElement(g_k,e);
////		FEMFunc fc = null;
////		if(g_c != null) fc = Utils.interpolateOnElement(g_c,e);
//
//		//所有自由度双循环
//		for(int i=1;i<=nDOFs;i++) {
//			DOF dofI = DOFs.at(i);
//			ScalarShapeFunction sfI = dofI.getSSF();
//			int nLocalRow = dofI.getLocalIndex();
//			int nGlobalRow = dofI.getGlobalIndex();
//			for(int j=1;j<=nDOFs;j++) {
//				DOF dofJ = DOFs.at(j);
//				int nLocalCol = dofJ.getLocalIndex();
//				int nGlobalCol = dofJ.getGlobalIndex();
//				//Integrand part of Weak Form on element e
//				FEMFunc integrand = null;
//				if(g_k == null) {
//					integrand = mapShape_x.get(nLocalRow).M(mapShape_x.get(nLocalCol))
//								.A(
//								mapShape_y.get(nLocalRow).M(mapShape_y.get(nLocalCol))
//								);
//				} else {
//					integrand = fk.M(
//									mapShape_x.get(nLocalRow).M(mapShape_x.get(nLocalCol))
//									.A(
//									mapShape_y.get(nLocalRow).M(mapShape_y.get(nLocalCol))
//									)
//								.A(
//									fc.M(dofI.getSSF().M(dofJ.getSSF())))
//								);
//				}
//				//Numerical integration on element e
//				double lhsVal = 0.0;
//				if(e.vertices().size() == 3) {
//					lhsVal = FOIntegrate.intOnTriangleRefElement(
//							integrand.M(e.getJacobin()),4
//							);
//				} else if (e.vertices().size() == 4) {
//					lhsVal = FOIntegrate.intOnRectangleRefElement(
//							integrand.M(e.getJacobin()),2 //TODO
//							);
//				}
//				globalStiff.add(nGlobalRow, nGlobalCol, lhsVal);
//			}
//			//Load vector
//			FEMFunc ff = Utils.interpolateOnElement(g_f, e);
//			FEMFunc integrand = ff.M(sfI);
//			double rhsVal = 0.0;
//			if(e.vertices().size() == 3) {
//				rhsVal = FOIntegrate.intOnTriangleRefElement(
//						integrand.M(e.getJacobin()),4
//						);
//			} else if (e.vertices().size() == 4) {
//				rhsVal = FOIntegrate.intOnRectangleRefElement(
//						integrand.M(e.getJacobin()),2 //TODO
//						);
//			}
//			globalLoad.add(nGlobalRow, rhsVal);
//		}
//		
////		//Robin:  d*u + k*u_n= q (自然边界：d==k, q=0)
////		//if(g_d != null && e.isBorderElement()) {
////		if(e.isBorderElement()) {
////
////			ElementList beList = e.getBorderElements();
////			for(int n=1;n<=beList.size();n++) {
////				Element be = beList.at(n);
////				
////				FEMFunc fd = null;
////				if(g_d != null) fd = Utils.interpolateOnElement(g_d, be);
////				//Check node type
////				NodeType nodeType = be.getBorderNodeType();
////				if(nodeType == NodeType.Neumann || nodeType == NodeType.Robin) {
////					DOFList beDOFs = be.getAllDOFList(DOFOrder.NEFV);
////					int nBeDOF = beDOFs.size();
////					
////					//Update Jacobin on be
////					be.updateJacobinLinear1D();
////					
////					//形函数计算需要和单元关联
////					for(int i=1;i<=nBeDOF;i++) {
////						beDOFs.at(i).getSSF().assignElement(be);
////					}
////					
////					//所有自由度双循环
////					for(int i=1;i<=nBeDOF;i++) {
////						DOF dofI = beDOFs.at(i);
////						ScalarShapeFunction sfI = dofI.getSSF();
////						int nGlobalRow = dofI.getGlobalIndex();
////						if(g_d != null) {
////							for(int j=1;j<=nBeDOF;j++) {
////								DOF dofJ = beDOFs.at(j);
////								ScalarShapeFunction sfJ = dofJ.getSSF();
////								int nGlobalCol = dofJ.getGlobalIndex();
////								//Stiff matrix for border
////								FEMFunc borderIntegrand = fd.M(sfI.M(sfJ));
////								//Numerical integrate the border 'be' of element 'e'
////								double lhsBrVal = FOIntegrate.intOnLinearRefElement(
////										borderIntegrand.M(be.getJacobin()),5
////									);
////								globalStiff.add(nGlobalRow, nGlobalCol, lhsBrVal);
////							}
////						}
////						//Load vector for border
////						if(g_g != null) {
////							FEMFunc fq = Utils.interpolateOnElement(g_g, be);
////							FEMFunc borderIntegrand = fq.M(sfI);
////							double rhsBrVal = FOIntegrate.intOnLinearRefElement(
////									borderIntegrand.M(be.getJacobin()),5
////								);
////							globalLoad.add(nGlobalRow, rhsBrVal);
////						}
////					}
////				}
////			}
////		}
//	}
}
