package edu.uta.futureye.lib.weakform;

import java.util.HashMap;
import java.util.Map;

import edu.uta.futureye.core.DOF;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.intf.WeakForm;
import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.basic.FX;
import edu.uta.futureye.function.intf.MathFun;
import edu.uta.futureye.function.intf.ScalarShapeFunction;
import edu.uta.futureye.function.intf.VectorShapeFunction;
import edu.uta.futureye.util.FutureyeException;
import edu.uta.futureye.util.Utils;

public class WeakFormBuilder {
	protected DOF _trialDOF = null;
	protected DOF _testDOF = null;
	protected ScalarShapeFunction _u = null;
	protected ScalarShapeFunction _v = null;
	protected VectorShapeFunction _uu = null;
	protected VectorShapeFunction _vv = null;
	protected int _uDOFLocalIndex;
	protected int _vDOFLocalIndex;
	
	public static enum Type {
		LHS_Domain, //Left Hand Side, integration on domain
		LHS_Border, //Left Hand Side, integration on Border
		RHS_Domain, //Right Hand Side, integration on domain
		RHS_Border  //Right Hand Side, integration on Border
		};
	
	Map<String, MathFun> param = new HashMap<String, MathFun>();
	
	public WeakFormBuilder(MathFun ...f) {
		for(int i=0;i<f.length;i++)
			this.addParamters(f[i]);
	}
	
	public WeakFormBuilder addParamters(MathFun f) {
		if(f.getFunName() == null)
			throw new FutureyeException("Please specify function name!");
		param.put(f.getFunName(), f);
		return this;
	}
	
	public WeakFormBuilder addParam(String fName, MathFun f) {
		param.put(fName, f);
		return this;
	}
	
	public MathFun getParam(String fName, Element e) {
		MathFun f = param.get(fName);
		if(f == null)
			throw new FutureyeException("Can NOT find the parameter "+fName+"!");
		MathFun ff = Utils.interpolateOnElement(f,e);
		return ff;
	}
	
	public ScalarShapeFunction getScalarTrial() {
		return _u;
	}
	
	public ScalarShapeFunction getScalarTest() {
		return _v;
	}
	
	public VectorShapeFunction getVectorTrial() {
		return _uu;
	}
	
	public VectorShapeFunction getVectorTest() {
		return _vv;
	}
	
	public int getTrialDOFLocalIndex() {
		return _uDOFLocalIndex;
	}
	
	public int getTestDOFLocalIndex() {
		return _vDOFLocalIndex;
	}
	public DOF getTrialDOF() {
		return _trialDOF;
	}
	
	public DOF getTestDOF() {
		return _testDOF;
	}
	
	/**
	 * Implement yourself weak form here
	 * <p><blockquote><pre>
	 * //the following example equals to class WeakFormLaplace2D
	 * ScalarShapeFunction u = getScalarTrial();
	 * ScalarShapeFunction v = getScalarTest();
	 * Function fk = param(e,"k");
	 * Function fc = param(e,"c");
	 * Function fd = param(e,"d");
	 * Function ff = param(e,"f");
	 * Function fq = param(e,"q");
	 * 
	 * switch(type) {
	 * case LHS_Domain:
	 * 		return fk.M(
	 * 				u._d("x").M(v._d("x")) .A (u._d("y").M(v._d("y")))
	 * 			).A(
	 * 				fc.M(u.M(v))
	 * 			);
	 * 	case LHS_Border:
	 * 		return fd.M(u.M(v));
	 * 	case RHS_Domain:
	 * 		return ff.M(v);
	 * 	case RHS_Border:
	 * 		return fq.M(v);
	 * 	default:
	 * 		return null;
	 * }
	 * </blockquote></pre>
	 * @param e
	 * @param type
	 * @return
	 */
	public MathFun makeExpression(Element e, Type type) {
		throw new FutureyeException("Override this function to implement your Weak Form!");
	}
	
	public WeakForm getScalarWeakForm() {
		WeakForm wf = new AbstractScalarWeakForm() {
			@Override
			public MathFun leftHandSide(Element e, ItemType itemType) {
				if(itemType==ItemType.Domain)  {
					//Integrand part of Weak Form on element e
					MathFun integrand = makeExpression(e,Type.LHS_Domain);
					return integrand;
				} else if(itemType==ItemType.Border) {//Neumann border integration on LHS
					MathFun borderIntegrand = makeExpression(e,Type.LHS_Border);
					return borderIntegrand;
				}
				return null;
			}

			@Override
			public MathFun rightHandSide(Element e, ItemType itemType) {
				if(itemType==ItemType.Domain)  {
					MathFun integrand = makeExpression(e,Type.RHS_Domain);
					return integrand;
				} else if(itemType==ItemType.Border) {
					MathFun borderIntegrand = makeExpression(e,Type.RHS_Border);
					return borderIntegrand;
				}
				return null;	
			}
			
			@Override
			public void setDOF(DOF trialDOF, DOF testDOF) {
				super.setDOF(trialDOF, testDOF);
				if(trialDOF != null) {
					_trialDOF = trialDOF;
					_u = trialDOF.getSSF();
					_uDOFLocalIndex = trialDOF.getLocalIndex();
				} 
				if(testDOF != null) {
					_testDOF = testDOF;
					_v = testDOF.getSSF();
					_vDOFLocalIndex = testDOF.getLocalIndex();
				}
			}
		};
		return wf;
	}
	
	public WeakForm getVectorWeakForm() {
		WeakForm wf = new AbstractVectorWeakForm() {
			@Override
			public MathFun leftHandSide(Element e, ItemType itemType) {
				if(itemType==ItemType.Domain)  {
					//Integrand part of Weak Form on element e
					MathFun integrand = makeExpression(e,Type.LHS_Domain);
					return integrand;
				} else if(itemType==ItemType.Border) {//Neumann border integration on LHS
					MathFun borderIntegrand = makeExpression(e,Type.LHS_Border);
					return borderIntegrand;
				}
				return null;
			}

			@Override
			public MathFun rightHandSide(Element e, ItemType itemType) {
				if(itemType==ItemType.Domain)  {
					MathFun integrand = makeExpression(e,Type.RHS_Domain);
					return integrand;
				} else if(itemType==ItemType.Border) {
					MathFun borderIntegrand = makeExpression(e,Type.RHS_Border);
					return borderIntegrand;
				}
				return null;	
			}		
			
			@Override
			public void setDOF(DOF trialDOF, DOF testDOF) {
				super.setDOF(trialDOF, testDOF);
				if(trialDOF != null) {
					_trialDOF = trialDOF;
					_uu = trialDOF.getVSF();
					_uDOFLocalIndex = trialDOF.getLocalIndex();
				} 
				if(testDOF != null) {
					_testDOF = testDOF;
					_vv = testDOF.getVSF();
					_vDOFLocalIndex = testDOF.getLocalIndex();
				}
			}
		};
		return wf;
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		WeakFormBuilder wfb = new WeakFormBuilder() {
			@Override
			public MathFun makeExpression(Element e, Type type) {
				//example: equals to class WeakFormLaplace2D
				ScalarShapeFunction u = getScalarTrial();
				ScalarShapeFunction v = getScalarTest();
				MathFun fk = getParam("k",e);
				MathFun fc = getParam("c",e);
				MathFun fd = getParam("d",e);
				MathFun ff = getParam("f",e);
				MathFun fq = getParam("q",e);
				
				switch(type) {
					case LHS_Domain:
						return fk.M(
								u._d("x").M(v._d("x")) .A (u._d("y").M(v._d("y")))
							).A(
								fc.M(u.M(v))
							);
					case LHS_Border:
						return fd.M(u.M(v));
					case RHS_Domain:
						return ff.M(v);
					case RHS_Border:
						return fq.M(v);
					default:
							return null;
				}
			}			
		};
		wfb.addParam("k", FC.C1);
		wfb.addParam("c", FC.C0);
		wfb.addParam("d", FC.C1);
		wfb.addParam("f", FX.fx.M(FX.fx));
		wfb.addParam("q", FC.C0);
		WeakForm wf = wfb.getScalarWeakForm();
		System.out.println(wf.getTrialDOF());

	}

}
