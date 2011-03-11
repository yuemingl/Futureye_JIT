package edu.uta.futureye.function.basic;

import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.core.Edge;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Face;
import edu.uta.futureye.core.geometry.GeoEntity;
import edu.uta.futureye.function.AbstractFunction;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.intf.ElementDependentFunction;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.function.operator.FOBasic;
import edu.uta.futureye.function.operator.FOVector;
import edu.uta.futureye.util.FutureyeException;

/**
 * u关于n的方向导数：
 * \frac{ \partial{u} }{ \partial{\mathbf{n}} }
 * 
 * @author liuyueming
 *
 */
public class DuDn extends AbstractFunction implements ElementDependentFunction {
	protected Element e = null;
	protected Function u = null;
	protected Function u_x = null;
	protected Function u_y = null;
	protected Function u_z = null;
	protected Vector norm = null;
	
	public DuDn(Function u) {
		this.u = u;
	}
	
	public DuDn(Function u_x, Function u_y, Function u_z) {
		this.u_x = u_x;
		this.u_y = u_y;
		this.u_z = u_z;
	}
	
	@Override
	public void setElement(Element e) {
		this.e = e;
		GeoEntity ge = e.getGeoEntity();
		if(ge instanceof Edge) {
			norm = ((Edge)ge).getNormVector();
		} else if(ge instanceof Face) {
			norm = ((Face)ge).getNormVector();
		} else {
			FutureyeException ex = new FutureyeException("Unsuported element type");
			ex.printStackTrace();
			System.exit(-1);
		}
	}

	@Override
	public double value(Variable v) {
		Function rlt = null;
		this.setElement(v.getElement());
		if(u != null) {
			rlt = FOVector.Grad(u).dot(norm);
		} else if(this.norm.getDim() == 2) {
			rlt = FOBasic.Plus(
					FOBasic.Mult(u_x, new FC(norm.get(1))),
					FOBasic.Mult(u_y, new FC(norm.get(2)))
					);
		} else if(this.norm.getDim() == 3) {
			rlt = FOBasic.PlusAll(
					FOBasic.Mult(u_x, new FC(norm.get(1))),
					FOBasic.Mult(u_y, new FC(norm.get(2))),
					FOBasic.Mult(u_z, new FC(norm.get(3)))
					);
		} else {
			FutureyeException ex = new FutureyeException("Error");
			ex.printStackTrace();
			System.exit(-1);
		}
		return rlt.value(v);
	}
	
	public String toString() {
		return "DuDn";
	}
}
