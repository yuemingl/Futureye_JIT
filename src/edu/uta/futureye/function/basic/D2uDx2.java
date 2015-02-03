package edu.uta.futureye.function.basic;

import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.function.AbstractMathFun;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.intf.ElementDependentFunction;
import edu.uta.futureye.function.intf.MathFun;
import edu.uta.futureye.util.FutureyeException;

/**
 * Compute derivative u_xx of a vector u that defined on a mesh
 * There are two types of return values of u_xx
 * 
 * <p>
 * 1. Continuous on the whole mesh
 * 
 * <p><blockquote><pre>
 * //Assume u is a vector defined on a mesh
 * Function fu = new Vector2Function(u,"x","y"); //Convert to a function class
 * D2uDx2 d2udx2 = new D2uDx2(mesh, fu, "x"); //Define second derivative d2udx2 of fu
 * //Evaluate d2udx2 at {@code node}
 * Variable var = new Variable.createFrom(u, node, node.globalIndex);
 * System.out.println(d2udx2.value(v));
 * </pre></blockquote>
 * <p>
 * 2. Continuous on each element e, but may be discontinuous on the whole mesh.
 *    This feature can be used in assembling process for equations with coefficient
 *    that have derivative forms in it.
 *    
 * <p><blockquote><pre>
 * //Assume u is a vector defined on a mesh
 * Function fu = new Vector2Function(u,"x","y"); //Convert to a function class
 * D2uDx2 d2udx2 = new D2uDx2(mesh, fu, "x"); //Define derivative dudx of fu
 * //Evaluate d2udx2 at {@code node} which is specified that defined on element e
 * Variable var = new Variable.createFrom(u, node, node.globalIndex);
 * v.setElement(e);
 * System.out.println(d2udx2.value(v));
 * </pre></blockquote>
 * 
 * <p> 
 * NOTE:
 * <p> 
 *   This class now supports the case of 2D bilinear element only.
 * 
 * @see D2uDx2
 * @author liuyueming
 *
 */
public class D2uDx2 extends AbstractMathFun implements ElementDependentFunction {
	protected Element e = null;
	protected Mesh mesh = null;
	protected Vector2Function u = null;
	protected String x = null;
	protected MathFun fdu2 = null;
	
	/**
	 * 
	 * @param mesh
	 * @param u: u(x,y) or u(x,y,z)
	 * @param x
	 */
	public D2uDx2(Mesh mesh, Vector2Function u, String x) {
		this.mesh = mesh;
		this.u = u;
		this.x = x;
		this.setVarNames(u.getVarNames());
	}
	
	@Override
	public void setElement(Element e) {
		this.e = e;
	}

	@Override
	public double apply(Variable v) {
		throw new FutureyeException("unsupported error!");
		/*
		Element newEle = v.getElement();
		if(newEle != null) {
			Element ve = v.getElement();
			//3D?
			if(ve.nodes.size()==2 && ve.getGeoEntity() instanceof Edge) {
				//find the cell that includes this edge
				Node n1 = ve.nodes.at(1);
				Node n2 = ve.nodes.at(2);
				int findIndex = -1;
				for(int i=1;i<=n1.belongToElements.size();i++) {
					for(int j=1;j<=n2.belongToElements.size();j++) {
						int idx1 = n1.belongToElements.at(i).globalIndex;
						int idx2 = n2.belongToElements.at(j).globalIndex;
						if(idx1 == idx2) {
							findIndex = idx1;
							break;
						}
						if(findIndex > 0) break;
					}
				}
				this.setElement(mesh.getElementList().at(findIndex));
			}
			
		
			int N = e.nodes.size();
			double[] f = new double[N];
			for(int i=1;i<=N;i++) {
				Node node = e.nodes.at(i);
				Variable var = Variable.createFrom(u, node, node.globalIndex);
				f[i-1] = u.value(var);
			}
			double[] a = Utils.computeBilinearFunctionCoef(e.nodes.toArray(new Point[0]), f);
			//d2(a1 + a2*x + a3*y + a4*x*y)/dx2
			//d2(a1 + a2*x + a3*y + a4*x*y)/dy2
//???
			
		} else {
			if(this.fdu2 == null) {
				Vector du2 = Tools.computeDerivativeFast(mesh, u.u, x);
				du2 = Tools.computeDerivativeFast(mesh, du2, x);
				fdu2 = new Vector2Function(du2);
			}
			return fdu2.value(v);
		}
		*/
	}
	
	public String toString() {
		return "D2uDx2";
	}
}
