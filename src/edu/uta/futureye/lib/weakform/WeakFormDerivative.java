package edu.uta.futureye.lib.weakform;


import edu.uta.futureye.core.DOF;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.basic.Vector2Function;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.function.intf.ScalarShapeFunction;
import edu.uta.futureye.util.container.DOFList;

/**
 * 用有限元方法求导数
 * Solve: (w, v) = (U_x, v)
 * where w is unknown
 *   U_x is the piecewise derivative on the mesh
 *   w is an approximation of U_x
 *   
 * @author liuyueming
 */
public class WeakFormDerivative extends AbstractScalarWeakForm {
	protected Vector2Function g_U = null;
	protected String varName; // "x" or "y"

	public WeakFormDerivative(String varName) {
		this.varName = varName;
	}
	
	@Override
	public Function leftHandSide(Element e, ItemType itemType) {
		if(itemType==ItemType.Domain)  {
			//Integrand part of Weak Form on element e
			Function integrand = u.M(v);
			return integrand;
		}
		return null;
	}

	@Override
	public Function rightHandSide(Element e, ItemType itemType) {
		if(itemType==ItemType.Domain)  {
			Function rlt = new FC(0.0);
			int nNode = e.nodes.size();
			for(int i=1;i<=nNode;i++) {
				DOFList dofListI = e.getNodeDOFList(i);
				for(int k=1;k<=dofListI.size();k++) {
					DOF dofI = dofListI.at(k);
					Variable var = Variable.createFrom(g_U, (Node)dofI.getOwner(), dofI.getGlobalIndex());
					Function PValue = new FC(g_U.value(var));
					ScalarShapeFunction shape = dofI.getSSF();
					//以前版本需要调用shapeFun.asignElement(e)，现在版本不需要调用了
					rlt = rlt.A(PValue.M(shape._d(varName)));
				}
			}
			
			Function integrand = rlt.M(v);
			return integrand;
		}
		return null;
	}

	public void setParam(Vector2Function U) {
		this.g_U = U;
	}
}
