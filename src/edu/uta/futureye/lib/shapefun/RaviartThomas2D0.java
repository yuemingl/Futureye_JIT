package edu.uta.futureye.lib.shapefun;

import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import edu.uta.futureye.algebra.SpaceVector;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.core.CoordinateTransform;
import edu.uta.futureye.core.Edge;
import edu.uta.futureye.core.EdgeLocal;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.function.AbstractFunction;
import edu.uta.futureye.function.AbstractVectorFunction;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.FXY;
import edu.uta.futureye.function.basic.SpaceVectorFunction;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.function.intf.ScalarShapeFunction;
import edu.uta.futureye.function.intf.VectorFunction;
import edu.uta.futureye.function.intf.VectorShapeFunction;
import edu.uta.futureye.function.operator.FOVector;
import edu.uta.futureye.util.FutureyeException;
import edu.uta.futureye.util.container.ObjList;

/**
 * Raviart-Thomas 2D0 triangle element
 * 
 * for j=1,2,3
 * \psi_{E_{j}} = \sigma_{j} \frac{|E_{j}|}{2|T|} (\mathbf{x}-P_{j})
 *
 * where
 *   \mathbf{x} \in T
 *   E_{j} be the edge of a triangle T opposite to its vertex P_{j}
 *   sigma_{j} = \nu_{j} \dot \nu_{E_{j}}
 *   \nu_{j} the unit normal of T along E_{j} 
 *   \nu_{E_{j}} the unit normal vector of E_{j} chosen with a global fixed orientation
 * 
 * Property:
 *   \psi_{E_{j}} \dot \nu_{j} = 1, along E_{j}
 *                             = 0, along other edge
 * Ref.
 * C. Bahriawati, Three matlab implementations of the lowest-order Raviart-Thomas 
 * mfem with a posteriori error control, Computational Methods in Applied Mathematics, 
 * Vol.5(2005), No.4, pp.333-361.

 * e.g.
 * 
 *    P3
 *    | \ 
 *  E2|  \ E1
 *    |   \
 *    |    \
 *   P1-----P2 
 *       E3
 *       
 * @author liuyueming
 *
 */
public class RaviartThomas2D0 extends AbstractVectorFunction implements VectorShapeFunction {
	int funIndex = 0;
	
	private VectorFunction funCompose = null;
	private VectorFunction funOuter = null;
	private ObjList<String> innerVarNames = null;

	public RaviartThomas2D0(int funID) {
		funIndex = funID - 1;
		varNames.add("x");
		varNames.add("y");
		innerVarNames = new ObjList<String>("x","y");
	}
	
	
	@Override
	public void asignElement(final Element e) {
		//Space Vector Function
		Function fx = new FXY(varNames,1,0);
		Function fy = new FXY(varNames,0,1);
		SpaceVectorFunction svf = new SpaceVectorFunction(fx,fy);
	
		EdgeLocal lEdge = e.edges().at(funIndex+1);
		Edge gEdge = lEdge.getGlobalEdge();
		double sigma = lEdge.getNormVector().dot(gEdge.getNormVector());
		//|T| = area
		double area = e.getElementArea();
		//|E_{j}| = edgeLength
		double edgeLength = gEdge.getEdgeLength();
		double coef = sigma * edgeLength / (2.0 * area);
		
		int[] nodeIndex = {3, 1, 2};
		SpaceVector v = new SpaceVector(e.nodes.at(nodeIndex[funIndex]).coords());
		this.funOuter = FOVector.ScalarProduct(coef, 
				FOVector.Minus(svf, new SpaceVectorFunction(v)));
		
		//复合函数
		Map<String, Function> fInners = new HashMap<String, Function>();
		List<String> varNamesInner = new LinkedList<String>();
		varNamesInner.add("r");
		varNamesInner.add("s");

		for(final String varName : varNames) {
			fInners.put(varName, new AbstractFunction(varNamesInner) {
				
				protected CoordinateTransform trans = new CoordinateTransform(2);
				
				public Function _d(String varName) {
					return null;
				}
				@Override
				public double value(Variable v) {
					//根据不同的varName给出不同的表达式
					//x = x(r,s,t)
					//y = y(r,s,t)
					//where t = 1 - r - s
					List<Function> transFun = trans.getTransformFunction(
							trans.getTransformLinear2DShapeFunction(e));
					if(varName.equals("x")) {
						return transFun.get(0).value(v);
					} else if(varName.equals("y")) {
						return transFun.get(1).value(v);
					} else {
						Exception e = new FutureyeException("Error!");
						e.printStackTrace();
					}
					return 0.0;
				}
			});
		}
		
		//使用复合函数构造形函数
		//funOuter.setVarNames(varNames); //!!!
		funCompose = FOVector.Compose(funOuter, fInners);
	}

	@Override
	public ScalarShapeFunction restrictTo(int funIndex) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Function dot(VectorFunction b) {
		return this.funCompose.dot(b);
	}

	@Override
	public Function get(int index) {
		return this.funCompose.get(index);
	}

	@Override
	public void set(int index, Function value) {
		throw new UnsupportedOperationException();
	}
	
	@Override
	public int getDim() {
		return this.funCompose.getDim();
	}

	@Override
	public void print() {
		// TODO Auto-generated method stub
	}

	@Override
	public Vector value(Variable v) {
		return this.funCompose.value(v);
	}

	@Override
	public ObjList<String> innerVarNames() {
		return innerVarNames;
	}
}
