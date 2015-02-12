package edu.uta.futureye.core.intf;


import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.core.DOF;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.function.FEMFunc;
import edu.uta.futureye.function.intf.MathFun;

public interface WeakForm {
	static enum ItemType {Domain, Border};
	
	//--- Common approach providing weak form interfaces to assembler-----

//	void setShapeFunction(ShapeFunction trial, int trialDofLocalIndex,
//			ShapeFunction test, int testDofLocalIndex);
	
	/**
	 * Set DOF objects to the weak form. These objects contain shape functions,
	 * local and global index of DOFs and geometry information of the geometry
	 * objects that possess the corresponding DOF objects
	 * <p>
	 * 直接传入DOF对象。对于一般的弱形式，只需要形函数就够了，
	 * 但是对于需要特殊处理的弱形式（例如：迎风方法需要知道自由度
	 * 对应的几何实体信息），仅提供形函数是不够的，因此该接口
	 * 允许将DOF对象直接传递给弱形式
	 */
	void setDOF(DOF trialDOF, DOF testDOF);
	
	/**
	 * Get degree of freedom(DOF) which contains trial shape function
	 * 
	 * @return
	 */
	DOF getTrialDOF();
	
	/**
	 * Get degree of freedom(DOF) which contains test shape function
	 * 
	 * @return
	 */
	DOF getTestDOF();
	
	/**
	 * Left hand side of the weak form
	 * 
	 * @param e 
	 * @param itemType
	 * @return
	 */
	FEMFunc leftHandSide(Element e, ItemType itemType);
	
	/**
	 * Right hand side of the weak form
	 * @param e
	 * @param itemType
	 * @return
	 */
	FEMFunc rightHandSide(Element e, ItemType itemType);
	
	/**
	 * Provide a pre-process function before calling leftHandSide(...) and rightHandSide(...) if necessary
	 * @param e
	 */
	void preProcess(Element e);

	//----------------------------------------------------------
	
	//--- Fast approach providing weak form interface to assembler-----
	/**
	 * Assemble element <code>e</code> here, instead of providing left hand side
	 * and right hand side.
	 * 
	 * @param e
	 * @param globalStiff (I/O): Global stiff matrix 
	 * @param globalLoad (I/O): Global load vector
	 *   
	 */
	void assembleElement(Element e, 
			Matrix globalStiff, Vector globalLoad);
	//--------------------------------------------------------
	
	/**
	 * Integrate on element <code>e</code>
	 * 
	 * @param e
	 * @param fun: LHS or RHS
	 * @return
	 */
	double integrate(Element e, MathFun fun);
	
	/**
	 * This interface has NO meaning for scalar valued problems.
	 * For vector valued problems, it is used to indicate if two components of
	 * vector function are coupled or independent variables. This information
	 * can be used to simplify the assembling process. For example, for 2D Stokes 
	 * problem (u v p): u and v are independent while u and p or v and p are coupled. 
	 * 
	 * @param nComponent1
	 * @param nComponent2
	 * @return
	 */
	boolean isVVFComponentCoupled(int nComponent1, int nComponent2);
	
	
}
