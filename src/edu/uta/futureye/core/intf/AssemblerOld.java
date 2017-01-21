package edu.uta.futureye.core.intf;

import edu.uta.futureye.algebra.intf.SparseMatrix;
import edu.uta.futureye.algebra.intf.SparseVector;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.function.intf.VectorMathFunc;

/**
 * Assembler Interface
 * 
 * 整体合成接口，最初在求解Laplace问题时（标量值函数）并没有抽象为接口。
 * 但是，当处理向量值函数时，可能会引入分块矩阵，这时整体合成的过程与以前
 * 有很大区别，需要以另外的方式实现。但是这两类问题具有相同的基本操作，
 * 因此将“整体合成”过程的基本操作抽象为接口。不同问题的合成器需要实现该
 * 通用接口，但实现细节由合成器自己完成。
 * 
 * 
 * @author liuyueming
 *
 */
public interface AssemblerOld {
	/**
	 * Assemble element's contributions to global stiffness matrix and global load vector
	 * <p> 
	 * 执行整体合成操作，将单元贡献合成到全局刚度矩阵和全局载荷向量中。
	 */
	void assemble();

	/**
	 * Return assembled global stiffness matrix
	 * @return
	 */
	SparseMatrix getStiffnessMatrix();
	
	/**
	 * Return assembled global load vector
	 * @return
	 */
	SparseVector getLoadVector();
	
	/**
	 * Impose Dirichlet boundary condition constraints for
	 * scalar valued problems. This function will affect
	 * global stiffness matrix and global load vector
	 * 
	 * @param diri
	 */
	void imposeDirichletCondition(MathFunc diri);
	
	/**
	 * Impose Dirichlet boundary condition constraints for 
	 * vector valued problems. This function will affect
	 * global stiffness matrix and global load vector
	 * 
	 * @param diri
	 */
	void imposeDirichletCondition(VectorMathFunc diri);
}
