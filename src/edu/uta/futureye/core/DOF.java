package edu.uta.futureye.core;

import edu.uta.futureye.function.intf.ShapeFunction;

/**
 * Degree Of Freedom
 * @author liuyueming
 *
 */
public class DOF {
	protected int localIndex;
	protected int globalIndex;
	protected ShapeFunction shapeFunction;
	protected Node owner;
	
	public DOF(int localIndex, int globalIndex, ShapeFunction shape) {
		this.localIndex = localIndex;
		this.globalIndex = globalIndex;
		this.shapeFunction = shape;
	}
	
	public int getLocalNumber() {
		return localIndex;
	}
	public void setLocalNumber(int localNumber) {
		this.localIndex = localNumber;
	}
	public int getGlobalNumber() {
		return globalIndex;
	}
	public void setGlobalNumber(int globalNumber) {
		this.globalIndex = globalNumber;
	}
	public ShapeFunction getShapeFunction() {
		return shapeFunction;
	}
	public void setShapeFunction(ShapeFunction shapeFunction) {
		this.shapeFunction = shapeFunction;
	}
	public Node getOwnerNode() {
		return owner;
	}
	public void setOwnerNode(Node node) {
		owner = node;
	}
}
