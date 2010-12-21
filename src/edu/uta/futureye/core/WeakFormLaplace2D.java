package edu.uta.futureye.core;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import edu.uta.futureye.algebra.Matrix;
import edu.uta.futureye.core.intf.WeakForm;
import edu.uta.futureye.function.DerivativeIndicator;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.function.intf.ShapeFunction;
import edu.uta.futureye.function.operator.FOBasic;
import edu.uta.futureye.function.operator.FOIntegrate;
import edu.uta.futureye.util.DOFList;
import edu.uta.futureye.util.ElementList;
import edu.uta.futureye.util.PairElementMatrix;
import edu.uta.futureye.util.Utils;

/**
 * Solve
 *   A(u, v) = (f, v)
 * where
 *   A(u, v) = (k*u_x, v_x) + (k*u_y, v_y) + (c*u, v)
 *   k = k(x,y)
 *   c = c(x,y)
 * 
 * @author liuyueming
 *
 */
public class WeakFormLaplace2D implements WeakForm {
	protected ShapeFunction u = null;
	protected ShapeFunction v = null;
	protected int uDOFLocalIndex;
	protected int vDOFLocalIndex;
	
	protected Function g_f = null;
	protected Function g_k = null;
	protected Function g_c = null;
	protected Function g_q = null;
	protected Function g_d = null;

	public void setF(Function f) {
		this.g_f = f;
	}
	
	//Robin: d*u + k*u_n = q
	public void setParam(Function k,Function c,Function q,Function d) {
		this.g_k = k;
		this.g_c = c;
		this.g_q = q;
		this.g_d = d;
	}

	@Override
	public Function leftHandSide(Element e, ItemType itemType) {
		u.asignElement(e);
		v.asignElement(e);

		if(itemType==ItemType.Domain)  {
			//Integrand part of Weak Form on element e
			DerivativeIndicator di_x = new DerivativeIndicator("x1");
			DerivativeIndicator di_y = new DerivativeIndicator("y1");
			Function integrand = null;
			if(g_k == null) {
				//System.out.println(u.derivative(di_x));
				integrand = FOBasic.Plus(
					FOBasic.Mult(u.derivative(di_x), v.derivative(di_x)),
					FOBasic.Mult(u.derivative(di_y), v.derivative(di_y))
					);
				//Variable v = new Variable();
				//v.set("r", 0.0);
				//v.set("s", 0.0);
				//System.out.println(u.derivative(di_x).value(v));
			} else {
				
				Function fk = Utils.interplateFunctionOnElement(g_k,e);
				Function fc = Utils.interplateFunctionOnElement(g_c,e);
				
				integrand = FOBasic.Plus(
							FOBasic.Mult(fk, FOBasic.Plus(
							FOBasic.Mult(u.derivative(di_x), v.derivative(di_x)),
							FOBasic.Mult(u.derivative(di_y), v.derivative(di_y))
						)),
							FOBasic.Mult(fc, FOBasic.Mult(u, v))
						);
			}

			//Numerical integration on element e
			Function integral = null;
			if(e.getVertexList().size() == 3) {
				integral = FOIntegrate.intOnTriangleRefElement(
						FOBasic.Mult(integrand, e.getJacobin()),5
						);
			} else if (e.getVertexList().size() == 4) {
				integral = FOIntegrate.intOnRectangleRefElement(
						FOBasic.Mult(integrand, e.getJacobin()),2 //TODO !!! 1
						);
			}
	
			return integral;
		}
		else if(itemType==ItemType.Border) {
			if(g_d != null) {
				Element be = e;
				Function fd = Utils.interplateFunctionOnElement(g_d, e);
				Function borderIntegrand = FOBasic.Mult(FOBasic.Mult(fd, u), v);
				//Numerical integration the border of element e
				//System.out.println(be.getJacobin().value(null));
				Function borderIntgral = FOIntegrate.intOnLinearRefElement(
						FOBasic.Mult(borderIntegrand, be.getJacobin()),5
					);
				return borderIntgral;
			}
		}
		return null;
	}

	@Override
	public Function rightHandSide(Element e, ItemType itemType) {
		v.asignElement(e);
		if(itemType==ItemType.Domain)  {
			Function ff = Utils.interplateFunctionOnElement(g_f, e);
			Function integrand = FOBasic.Mult(ff,v);
			
			Function integral = null;
			if(e.getVertexList().size() == 3) {
				integral = FOIntegrate.intOnTriangleRefElement(
						FOBasic.Mult(integrand, e.getJacobin()),5
					);
			} else if (e.getVertexList().size() == 4) {
				integral = FOIntegrate.intOnRectangleRefElement(
						FOBasic.Mult(integrand, e.getJacobin()),2 //TODO !!! 1
						);
			}
			//System.out.println(e.nodes.size()+":  "+integral.value(null));
			return integral;
		} else if(itemType==ItemType.Border) {//add 2010-9-28
			Element be = e;
			Function fq = Utils.interplateFunctionOnElement(g_q, e);
			Function borderIntegrand = FOBasic.Mult(fq, v);
			//Numerical integration the border of element e
			//System.out.println(be.getJacobin().value(null));
		
			Function borderIntgral = FOIntegrate.intOnLinearRefElement(
					FOBasic.Mult(borderIntegrand, be.getJacobin()),5
				);
			//System.out.println(borderIntgral.value(null));
			return borderIntgral;
		} else
			return null;	
	}

	/**
	 * N = N(r,s) = N( r(x,y), s(x,y) )
	 */
	@Override
	public void setShapeFunction(ShapeFunction trial, int trialDofLocalIndex,
			ShapeFunction test, int testDofLocalIndex) {
		u = trial;
		v = test;
		this.uDOFLocalIndex = trialDofLocalIndex;
		this.vDOFLocalIndex = testDofLocalIndex;
	}

	/**
	 * 返回的两个矩阵怎么处理?
	 */
	@Override
//	public List<PairElementMatrix> associateElement(Element e) {
//		return null;
//	}
	public List<PairElementMatrix> associateElement(Element e) {
		
		
		List<PairElementMatrix> retList= new ArrayList<PairElementMatrix>();
		
		e.updateJacobinLinear2D();
		int nNode = e.nodes.size();
		
		Map<Integer, Function> mapShape_x = new LinkedHashMap<Integer, Function>();
		Map<Integer, Function> mapShape_y = new LinkedHashMap<Integer, Function>();
		DerivativeIndicator di_x = new DerivativeIndicator("x1");
		DerivativeIndicator di_y = new DerivativeIndicator("y1");
		
		//associate element to shape function
		for(int i=1;i<=nNode;i++) {
			DOFList dofListI = e.getDOFList(i);
			int nDOF = dofListI.size();
			for(int k=1;k<=nDOF;k++) {
				DOF dof = dofListI.at(k);
				dof.shapeFunction.asignElement(e);
				mapShape_x.put(dof.localIndex, dof.shapeFunction.derivative(di_x));
				mapShape_y.put(dof.localIndex, dof.shapeFunction.derivative(di_y));
			}
		}
		
		int dofDim = e.getTotalNumberOfDOF();
		PairElementMatrix pem = 
			new PairElementMatrix(e,new Matrix(dofDim,dofDim));

		Function fk = null;
		if(g_k != null) fk = Utils.interplateFunctionOnElement(g_k,e);
		Function fc = null;
		if(g_c != null) fc = Utils.interplateFunctionOnElement(g_c,e);
		Function fd = null;
		if(g_d != null) fd = Utils.interplateFunctionOnElement(g_d, e);

		//结点双循环
		for(int i=1;i<=nNode;i++) {
			for(int j=1;j<=nNode;j++) {
				DOFList dofListI = e.getDOFList(i);
				DOFList dofListJ = e.getDOFList(j);
				
				//结点上自由度双循环
				//不同结点上的自由度个数不一定相同(e.g. adaptive finite element)
				int nDOF1 = dofListI.size();
				int nDOF2 = dofListJ.size();
				for(int k1=1;k1<=nDOF1;k1++) {
					DOF dofI = dofListI.at(k1);
					for(int k2=1;k2<=nDOF2;k2++) {
						DOF dofJ = dofListJ.at(k2);
						int lRow = dofI.localIndex;
						int lCol = dofJ.localIndex;
						
						//Integrand part of Weak Form on element e
						Function integrand = null;
						if(g_k == null) {
							integrand = FOBasic.Plus(
								FOBasic.Mult(mapShape_x.get(lRow), mapShape_x.get(lCol)),
								FOBasic.Mult(mapShape_y.get(lRow), mapShape_y.get(lCol))
								);
						} else {
							integrand = FOBasic.Plus(
										FOBasic.Mult(fk, FOBasic.Plus(
										FOBasic.Mult(mapShape_x.get(lRow), mapShape_x.get(lCol)),
										FOBasic.Mult(mapShape_y.get(lRow), mapShape_y.get(lCol))
									)),
										FOBasic.Mult(fc, FOBasic.Mult(dofI.shapeFunction, dofJ.shapeFunction))
									);
						}

						//Numerical integration on element e
						Function integral = null;
						if(e.getVertexList().size() == 3) {
							integral = FOIntegrate.intOnTriangleRefElement(
									FOBasic.Mult(integrand, e.getJacobin()),5
									);
						} else if (e.getVertexList().size() == 4) {
							integral = FOIntegrate.intOnRectangleRefElement(
									FOBasic.Mult(integrand, e.getJacobin()),2 //TODO !!! 1
									);
						}
						double rlt = integral.value(null);
						//System.out.println(rlt);
						pem.localMatrix.set(lRow, lCol, rlt);
					}
				}
			}
		}
		retList.add(pem);
		
		if(e.isBorderElement()) {
			ElementList beList = e.getSubElements();
			for(int ibeList=1;ibeList<=beList.size();ibeList++) {
				Element be = beList.at(ibeList);
				dofDim = be.getTotalNumberOfDOF();
				Edge edge = (Edge)be.stemGeoEntity;

				if(edge.getBorderNodeType() == NodeType.Robin) {
					be.updateJacobinLinear1D();
					nNode = be.nodes.size();
					//bugfix 12/27/2010 e->be
					PairElementMatrix pem2 = 
						new PairElementMatrix(be,new Matrix(dofDim,dofDim));
					for(int i=1;i<=nNode;i++) {
						for(int j=1;j<=nNode;j++) {
							DOFList dofListI = be.getDOFList(i);
							DOFList dofListJ = be.getDOFList(j);
							
							int nDOF1 = dofListI.size();
							int nDOF2 = dofListJ.size();
							for(int k1=1;k1<=nDOF1;k1++) {
								DOF dofI = dofListI.at(k1);
								for(int k2=1;k2<=nDOF2;k2++) {
									DOF dofJ = dofListJ.at(k2);

									int lRow = dofI.localIndex;
									int lCol = dofJ.localIndex;
									
									if(g_d != null) {
										Function borderIntegrand = FOBasic.Mult(
												FOBasic.Mult(fd, dofI.shapeFunction), dofJ.shapeFunction
												);
										//Numerical integrate the border 'be' of element 'e'
										//System.out.println(be.getJacobin().value(null));
										Function borderIntgral = FOIntegrate.intOnLinearRefElement(
												FOBasic.Mult(borderIntegrand, be.getJacobin()),5
											);
										double rlt = borderIntgral.value(null);
										//System.out.println(rlt);
										pem2.localMatrix.plusValue(lRow, lCol, rlt);
									}
								}
							}
						}
					}
					retList.add(pem2);
				}
			}
		}
		return retList;
	}
}
