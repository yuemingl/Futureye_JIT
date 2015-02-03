package edu.uta.futureye.function.basic;

import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.core.Edge;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Face;
import edu.uta.futureye.core.geometry.GeoEntity;
import edu.uta.futureye.function.AbstractMathFun;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.intf.ElementDependentFunction;
import edu.uta.futureye.function.intf.MathFun;
import edu.uta.futureye.function.operator.FMath;
import edu.uta.futureye.util.FutureyeException;

/**
 * u关于n的方向导数：
 * \frac{ \partial{u} }{ \partial{\mathbf{n}} }
 * 
 * @author liuyueming
 *
 */
public class DuDn extends AbstractMathFun implements ElementDependentFunction {
	protected Element e = null;
	protected MathFun u = null;
	protected MathFun u_x = null;
	protected MathFun u_y = null;
	protected MathFun u_z = null;
	protected Vector norm = null;
	
	public DuDn(MathFun u) {
		this.u = u;
	}
	
	public DuDn(MathFun u_x, MathFun u_y, MathFun u_z) {
		this.u_x = u_x;
		this.u_y = u_y;
		this.u_z = u_z;
	}
	
	@Override
	public void setElement(Element e) {
		this.e = e;
		//Compute outer normal vector 
		GeoEntity ge = e.getGeoEntity();
		if(ge instanceof Edge) {
			norm = ((Edge)ge).getNormVector();
		} else if(ge instanceof Face) {
			norm = ((Face)ge).getNormVector();
		} else {
			throw new FutureyeException("Unsuported element type");
		}
	}

	@Override
	public double apply(Variable v) {
		MathFun rlt = null;
		this.setElement(v.getElement());
		if(u != null) {
			//u is passed into constructor
			rlt = FMath.grad(u).dot(norm);
		} else if(this.norm.getDim() == 2) {
			//2D case
			rlt = u_x.M(new FC(norm.get(1)))
					.A(
				  u_y.M(new FC(norm.get(2)))
					);
		} else if(this.norm.getDim() == 3) {
			//3D case
			rlt = FMath.sum(
					u_x.M(new FC(norm.get(1))),
					u_y.M(new FC(norm.get(2))),
					u_z.M(new FC(norm.get(3)))
					);
		} else {
			throw new FutureyeException(
					"Error: u="+u+", this.norm.getDim()="+this.norm.getDim());
		}
		return rlt.apply(v);
	}
	
	public String toString() {
		return "DuDn";
	}
}
