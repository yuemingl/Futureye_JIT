package edu.uta.futureye.lib.weakform;


import edu.uta.futureye.core.DOF;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.basic.Vector2Function;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.function.intf.ScalarShapeFunction;
import edu.uta.futureye.function.operator.FOBasic;
import edu.uta.futureye.util.Utils;
import edu.uta.futureye.util.container.DOFList;

/**
 * 系数反问题
 * Solve: (U*u, v) = (f, v) - (k*grad(U),grad(v))
 * where u is unknown
 * U,f and k is known
 * 
 * @author liuyueming
 */
public class WeakFormL22D extends AbstractScalarWeakForm {
	protected Vector2Function g_U = null;
	protected Function g_f = null;
	protected Function g_k = null;


	@Override
	public Function leftHandSide(Element e, ItemType itemType) {
		if(itemType==ItemType.Domain)  {
			//Integrand part of Weak Form on element e
			Function integrand = null;
			Function fU = Utils.interplateFunctionOnElement(g_U, e);
			integrand = FOBasic.Mult(fU, FOBasic.Mult(u, v));
			return integrand;
		}
		return null;
	}

	@Override
	public Function rightHandSide(Element e, ItemType itemType) {
		if(itemType==ItemType.Domain)  {
			Function ff = Utils.interplateFunctionOnElement(g_f, e);
			Function fk = Utils.interplateFunctionOnElement(g_k, e);
			
			Function rlt_dx = new FC(0.0);
			Function rlt_dy = new FC(0.0);
			int nNode = e.nodes.size();
			for(int i=1;i<=nNode;i++) {
				DOFList dofListI = e.getNodeDOFList(i);
				for(int k=1;k<=dofListI.size();k++) {
					DOF dofI = dofListI.at(k);
					Variable var = Variable.createFrom(g_U, (Node)dofI.getOwner(), dofI.getGlobalIndex());
					Function PValue = new FC(g_U.value(var));
					ScalarShapeFunction shape = dofI.getSSF();
					//以前版本需要调用shapeFun.asignElement(e)，现在版本不需要调用了
					rlt_dx = FOBasic.Plus(rlt_dx, FOBasic.Mult(PValue, shape._d("x")));
					rlt_dy = FOBasic.Plus(rlt_dy, FOBasic.Mult(PValue, shape._d("y")));
				}
			}
			
			Function integrand = FOBasic.Minus(
				FOBasic.Mult(ff,v),
				FOBasic.Mult(fk,
						FOBasic.Plus(
								FOBasic.Mult(rlt_dx, v._d("x")),
								FOBasic.Mult(rlt_dy, v._d("y"))
			)));
			return integrand;
		}
		return null;
	}

	public void setF(Function f) {
		this.g_f = f;
	}
	
	public void setParam(Function k,Vector2Function U) {
		this.g_k = k;
		this.g_U = U;
	}

}
