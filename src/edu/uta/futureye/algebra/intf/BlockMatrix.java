/**
 * Copyright (c) 2010, nkliuyueming@gmail.com. All rights reserved.
 * 
 * 
 */
package edu.uta.futureye.algebra.intf;

import java.util.Map;

/**
 * <blockquote><pre>
 * General block matrix interface
 * 
 * e.g. A,B,C,D are all 3*3 matrix,
 * M = (A  B)
 *     (C  D)  
 * M.getBlock(2,1) == C
 * 
 * M as a general matrix,
 * M.get(1,4) == B(1,1)
 * </blockquote></pre>
 * 
 * @author liuyueming
 *
 */
public interface BlockMatrix<TMatrix extends Matrix> extends Matrix {
	
	public int getRowBlockDim();
	public int getColBlockDim();
	
	public TMatrix getBlock(int row, int col);
	public void setBlock(int row, int col,TMatrix m);
	
	public Map<Integer,Map<Integer,TMatrix>> getAllBlock();
	
}
