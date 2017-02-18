/**
 * Copyright (c) 2010, nkliuyueming@gmail.com. All rights reserved.
 * 
 * 
 */
package edu.uta.futureye.algebra.intf;

import java.util.Map;

/**
 * <blockquote><pre>
 * General block vector interface
 * 
 * e.g. A=(1 2 3)',B=(4 5)'
 * V = (A  B)'
 * V.getBlock(2) == B
 * 
 * V as a general vector,
 * V.get(4) == B(1)
 * </blockquote></pre>
 * 
 * @author liuyueming
 *
 */
public interface BlockVector<TVector extends Vector> extends Vector{
	
	public int getBlockDim();
	
	public TVector getBlock(int index);
	public void setBlock(int index,TVector v);
	
	public Map<Integer,TVector> getAllBlock();

}
