package edu.uta.futureye.function.basic;

import java.lang.reflect.Array;

import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.application.Tools;
import edu.uta.futureye.core.Edge;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Face;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.core.geometry.GeoEntity;
import edu.uta.futureye.core.geometry.Point;
import edu.uta.futureye.function.AbstractFunction;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.intf.ElementDependentFunction;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.function.operator.FMath;
import edu.uta.futureye.util.FutureyeException;
import edu.uta.futureye.util.Utils;

/**
 * Compute derivative u_x of a vector u that defined on a mesh
 * There are two types of return values of u_x
 * <p>
 * 1. Continuous on the whole mesh</br></br>
 * <code>
 * 	//Assume u is a vector define on a mesh</br>
 *  Function fu = new Vector2Function(u,"x","y"); //Convert to a function class</br>
 *	DuDx dudx = new DuDx(mesh, fu, "x"); //Define derivative dudx of fu</br>
 *  //Evaluate dudx at 'node'</br>
 *	Variable var = new Variable.createFrom(u, node, node.globalIndex);</br>
 *	System.out.println(dudx.value(v));</br>
 * </code></p>
 * <p>
 * 2. Continuous on each element e, but may be discontinuous on the whole mesh.</br>
 *    This feature can be used in assembling process for equations with coefficient
 *    that have derivative forms in it.</br></br>
 * <code>
 *  //Assume u is a vector define on a mesh</br>
 *  Function fu = new Vector2Function(u,"x","y"); //Convert to a function class</br>
 *	DuDx dudx = new DuDx(mesh, fu, "x"); //Define derivative dudx of fu</br>
 *  //Evaluate dudx at 'node' which is specified that defined on element e</br>
 *	Variable var = new Variable.createFrom(u, node, node.globalIndex);</br>
 *  v.setElement(e);</br>
 *	System.out.println(dudx.value(v));</br>
 * </code></p>
 * 
 * NOTE:
 *   This function supports 2D bilinear element case only now.
 * 
 * @author liuyueming
 *
 */
public class DuDx extends AbstractFunction implements ElementDependentFunction {
	protected Element e = null;
	protected Mesh mesh = null;
	protected Vector2Function u = null;
	protected String x = null;
	protected Function fdu2 = null;
	
	/**
	 * 
	 * @param mesh
	 * @param u: u(x,y) or u(x,y,z)
	 * @param x
	 */
	public DuDx(Mesh mesh, Vector2Function u, String x) {
		this.mesh = mesh;
		this.u = u;
		this.x = x;
		this.setVarNames(u.varNames());
		Array a;
	}
	
	@Override
	public void setElement(Element e) {
		this.e = e;
	}

	@Override
	public double value(Variable v) {
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
			//d(a1 + a2*x + a3*y + a4*x*y)/dx
			//d(a1 + a2*x + a3*y + a4*x*y)/dy
			Function du = null;
			if(x.equals("x")) {
				du = new FXY(0.0,a[3],a[1]);
			} else if(x.equals("y")) {
				du = new FXY(a[3],0.0,a[2]);
			} else {
				throw new FutureyeException("x(="+x+") should be 'x' or 'y'!");
			}
			Variable vv = new Variable();
			if(varNames==null || varNames.size()==0) {
				throw new FutureyeException("varNames should be specified in parameter 'u'.");
			}
			if(v.getIndex() != 0) {
				Node node = mesh.getNodeList().at(v.getIndex());
				for(int i=0;i<varNames.size();i++) {
					vv.set(varNames.get(i),node.coord(i+1));
				}
				return du.value(vv);
			} else {
				return du.value(v);	
			}
			
		} else {
			if(this.fdu2 == null) {
				Vector du2 = Tools.computeDerivativeFast(mesh, u.u, x);
				fdu2 = new Vector2Function(du2);
			}
			return fdu2.value(v);
		}
	}
	
	public String toString() {
		return "DuDn";
	}
}
