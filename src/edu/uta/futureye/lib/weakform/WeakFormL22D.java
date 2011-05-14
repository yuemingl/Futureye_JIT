package edu.uta.futureye.lib.weakform;


import edu.uta.futureye.core.DOF;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.basic.Vector2Function;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.function.intf.ScalarShapeFunction;
import edu.uta.futureye.util.Utils;
import edu.uta.futureye.util.container.DOFList;

/**
 * 系数反问题
 * Solve: (Uu, v) = (f, v) - (k\nabla{U},\nabla{v})
 *   where u is unknown
 *   U,f and k is known
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
			Function fU = Utils.interpolateFunctionOnElement(g_U, e);
			integrand = fU.M(u.M(v));
			return integrand;
		}
		return null;
	}

	@Override
	public Function rightHandSide(Element e, ItemType itemType) {
		if(itemType==ItemType.Domain)  {
			Function ff = Utils.interpolateFunctionOnElement(g_f, e);
			Function fk = Utils.interpolateFunctionOnElement(g_k, e);
			
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
					rlt_dx = rlt_dx.A(PValue.M(shape._d("x")));
					rlt_dy = rlt_dy.A(PValue.M(shape._d("y")));
				}
			}
			
			Function integrand = 
					ff.M(v)
					.S(
						fk.M(
							rlt_dx.M(v._d("x")).A(rlt_dy.M(v._d("y")))
						)
					);
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
