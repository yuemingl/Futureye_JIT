package edu.uta.futureye.application;

import edu.uta.futureye.algebra.intf.SparseMatrix;
import edu.uta.futureye.algebra.intf.SparseVector;

/**
 * Contain matrix and RHS of equation A*x=f
 * 
 * @author liuyueming
 *
 */
public class Equation {
	public SparseMatrix A;
	public SparseVector f;
}
