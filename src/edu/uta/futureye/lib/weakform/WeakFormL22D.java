package edu.uta.futureye.lib.weakform;


import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.core.geometry.Point;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.FXY;
import edu.uta.futureye.function.intf.MathFun;
import edu.uta.futureye.util.Utils;

/**
 * 系数反问题
 * Solve: (Uu, v) = (f, v) - (k\nabla{U},\nabla{v})
 *   where u is unknown
 *   U,f and k is known
 * 
 * @author liuyueming
 */
public class WeakFormL22D extends AbstractScalarWeakForm {
	protected MathFun g_U = null;
	protected MathFun g_Ux = null;
	protected MathFun g_Uy = null;
	protected MathFun g_f = null;
	protected MathFun g_k = null;


	@Override
	public MathFun leftHandSide(Element e, ItemType itemType) {
		if(itemType==ItemType.Domain)  {
			//Integrand part of Weak Form on element e
			MathFun integrand = null;
			MathFun fU = Utils.interpolateOnElement(g_U, e);
			integrand = fU.M(u.M(v));
			return integrand;
		}
		return null;
	}

	@Override
	public MathFun rightHandSide(Element e, ItemType itemType) {
		if(itemType==ItemType.Domain)  {
			MathFun ff = Utils.interpolateOnElement(g_f, e);
			MathFun fk = Utils.interpolateOnElement(g_k, e);
			
			MathFun integrand = null;
			if(g_Ux == null) {
//新方法1：计算导数				
				int N = e.nodes.size();
				double[] f = new double[N];
				for(int i=1;i<=N;i++) {
					Node node = e.nodes.at(i);
					Variable var = Variable.createFrom(g_U, node, node.globalIndex);
					f[i-1] = g_U.apply(var);
				}
				double[] a = Utils.computeBilinearFunctionCoef(e.nodes.toArray(new Point[0]), f);
				//d(a1 + a2*x + a3*y + a4*x*y)/dx
				MathFun dx = new FXY(0.0,a[3],a[1]);
				//d(a1 + a2*x + a3*y + a4*x*y)/dy
				MathFun dy = new FXY(a[3],0.0,a[2]);
				
				MathFun fUx = Utils.interpolateOnElement(dx, e);
				MathFun fUy = Utils.interpolateOnElement(dy, e);
				integrand = 
					ff.M(v)
					.S(
						fk.M(
							fUx.M(v._d("x")).A(fUy.M(v._d("y")))
						)
					);	
//新方法2：计算导数
				//利用类DuDx
				
				
//旧方法：计算导数				
//				Function rlt_dx = new FC(0.0);
//				Function rlt_dy = new FC(0.0);
//				int nNode = e.nodes.size();
//				for(int i=1;i<=nNode;i++) {
//					DOFList dofListI = e.getNodeDOFList(i);
//					for(int k=1;k<=dofListI.size();k++) {
//						DOF dofI = dofListI.at(k);
//						Variable var = Variable.createFrom(g_U, (Node)dofI.getOwner(), dofI.getGlobalIndex());
//						Function PValue = new FC(g_U.value(var));
//						ScalarShapeFunction shape = dofI.getSSF();
//						//以前版本需要调用shapeFun.asignElement(e)，现在版本不需要调用了
//						rlt_dx = rlt_dx.A(PValue.M(shape._d("x")));
//						rlt_dy = rlt_dy.A(PValue.M(shape._d("y")));
//					}
//				}
//				
//				integrand = 
//						ff.M(v)
//						.S(
//							fk.M(
//								rlt_dx.M(v._d("x")).A(rlt_dy.M(v._d("y")))
//							)
//						);
			} else {
				MathFun fUx = Utils.interpolateOnElement(g_Ux, e);
				MathFun fUy = Utils.interpolateOnElement(g_Uy, e);
				integrand = 
					ff.M(v)
					.S(
						fk.M(
							fUx.M(v._d("x")).A(fUy.M(v._d("y")))
						)
					);
				
			}
			return integrand;
		}
		return null;
	}

	public void setF(MathFun f) {
		this.g_f = f;
	}
	
	public void setParam(MathFun k,MathFun U) {
		this.g_k = k;
		this.g_U = U;
	}
	
	public void setParam(MathFun k, MathFun U, MathFun Ux, MathFun Uy) {
		this.g_k = k;
		this.g_U = U;
		this.g_Ux = Ux;
		this.g_Uy = Uy;
	}
}
