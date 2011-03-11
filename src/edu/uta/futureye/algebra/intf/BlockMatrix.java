package edu.uta.futureye.algebra.intf;

import java.util.Map;

/**
 * General block matrix interface
 * 
 * e.g. A,B,C,D are all 3*3 matrix,
 * M = (A  B)
 *     (C  D)  
 * M.getBlock(2,1) == C
 * 
 * M as a general matrix,
 * M.get(1,4) == B(1,1)
 * 
 * @author liuyueming
 *
 */
public interface BlockMatrix extends Matrix {
	
	public int getRowBlockDim();
	public int getColBlockDim();
	
	public Matrix getBlock(int row, int col);
	public void setBlock(int row, int col,Matrix m);
	
	public Map<Integer,Map<Integer,Matrix>> getAllBlock();
	
}
