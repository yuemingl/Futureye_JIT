package edu.uta.futureye.core;

import edu.uta.futureye.core.geometry.GeoEntity;
import edu.uta.futureye.function.intf.ScalarShapeFunction;
import edu.uta.futureye.function.intf.ShapeFunction;
import edu.uta.futureye.function.intf.VectorShapeFunction;

/**
 * Degree Of Freedom
 * 
 * U   = \sum_{j=1}^{N}(u_j*w), on whole domain
 * U_e = \sum_{i=1}^{n}(u_i*v), on element e
 * 
 * 
 * u_i*v => DOF object: <tt>dof</tt>
 * i     => <tt>dof.localIndex</tt>
 * v     => <tt>dof.shapeFun</tt>
 * j     => <tt>dof.globalIndex</tt>
 * w     => base function that constructed by combining
 *          shape functions on neighboring elements
 * 
 * @author liuyueming
 *
 */
public class DOF {
	protected int localIndex;
	protected int globalIndex;
	protected ShapeFunction shapeFun;
	protected GeoEntity owner;
	
	/**
	 * DOF belongs to component of vector valued function
	 * 	e.g. 2D Stokes Problem: (u v p)
	 *  nVVFComponent=1: DOF of u
	 *  nVVFComponent=2: DOF of v
	 *  nVVFComponent=3: DOF of p
	 */
	protected int nVVFComponent;
	
	/**
	 * Construct a Degree of Freedom(DOF) object
	 * 
	 * @param localIndex Local index of DOF
	 * @param globalIndex Global index of DOF
	 * @param shape Shape function associated with DOF
	 */
	public DOF(int localIndex, int globalIndex, ShapeFunction shape) {
		this.localIndex = localIndex;
		this.globalIndex = globalIndex;
		this.shapeFun = shape;
	}
	
	public int getLocalIndex() {
		return localIndex;
	}
	
	public void setLocalIndex(int localIndex) {
		this.localIndex = localIndex;
	}
	
	public int getGlobalIndex() {
		return globalIndex;
	}
	
	public void setGlobalIndex(int globalIndex) {
		this.globalIndex = globalIndex;
	}
	
	public int getVVFComponent() {
		return nVVFComponent;
	}
	public void setVVFComponent(int nComponent) {
		this.nVVFComponent = nComponent;
	}
	
	/**
	 * 返回形函数
	 * @return
	 */
	public ShapeFunction getSF() {
		return shapeFun;
	}
	
	/**
	 * 返回标量形函数
	 * @return
	 */
	public ScalarShapeFunction getSSF() {
		return (ScalarShapeFunction)shapeFun;
	}
	
	/**
	 * 返回向量形函数
	 * @return
	 */
	public VectorShapeFunction getVSF() {
		return (VectorShapeFunction)shapeFun;
	}

	
	public void setShapeFunction(ShapeFunction shapeFunction) {
		this.shapeFun = shapeFunction;
	}
	
	public GeoEntity getOwner() {
		return owner;
	}
	public void setOwner(GeoEntity enty) {
		owner = enty;
	}
	
	public String toString() {
		return String.format("{DOFL%d(G%d) SF:%s Owner:%s}", localIndex,
				globalIndex,
				shapeFun,
				owner);
	}
	
	public Node getNodeOwner() {
		return (Node)owner;
	}
	
	public EdgeLocal getLocalEdgeOwner() {
		return (EdgeLocal)owner;
	}
	
	public FaceLocal getLocalFaceOwner() {
		return (FaceLocal)owner;
	}
}
