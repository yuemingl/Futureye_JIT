package edu.uta.futureye.algebra.intf;

import java.util.Map;

/**
 * General block vector interface
 * 
 * e.g. A=(1 2 3)',B=(4 5)'
 * V = (A  B)'
 * V.getBlock(2) == B
 * 
 * V as a general vector,
 * V.get(4) == B(1)
 * 
 * @author liuyueming
 *
 */
public interface BlockVector extends Vector{
	
	public int getBlockDim();
	
	public Vector getBlock(int index);
	public void setBlock(int index,Vector v);
	
	public Map<Integer,Vector> getAllBlock();

}
