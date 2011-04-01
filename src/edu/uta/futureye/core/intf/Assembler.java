package edu.uta.futureye.core.intf;

import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.function.intf.VectorFunction;

/**
 * 整体合成接口，最初在求解Laplace问题时（标量值函数）并没有抽象为接口。
 * 但是，当处理向量值函数时，可能会引入分块矩阵，这时整体合成的过程与以前
 * 有很大区别，需要以另外的方式实现，因此将“整体合成”部分抽象为接口。
 * 
 * @author liuyueming
 *
 */
public interface Assembler {
	/**
	 * 执行整体合成操作
	 */
	void assemble();

	Matrix getStiffnessMatrix();
	
	Vector getLoadVector();
	
	void imposeDirichletCondition(Function diri);
	void imposeDirichletCondition(VectorFunction diri);
}
