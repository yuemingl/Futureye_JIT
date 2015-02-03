package edu.uta.futureye.function;

import edu.uta.futureye.function.intf.MathFun;
import edu.uta.futureye.function.intf.SimpleFun;
import edu.uta.futureye.function.intf.SimpleFun1;
import edu.uta.futureye.function.intf.SimpleFun2;
import edu.uta.futureye.function.intf.SimpleFun3;
import edu.uta.futureye.util.Constant;
import edu.uta.futureye.util.FutureyeException;

public class MathFunTODO extends AbstractMathFun {
	

	/**
	 * Construct a Function object. Lambda expression can be assigned to the functional interface SimpleFun. 
	 *
	 * 
	 * @param funDef: The definition of function itself
	 * @param funDerivatives: The first order derivatives of the function with respect to the free variables
	 */
	public MathFunTODO(SimpleFun funDef,
			SimpleFun ...funDerivatives) {
	}
	
	/**
	 * Construct a math function with one variable with undefined derivative
	 * 
	 * @param funDef
	 */
	public MathFunTODO(SimpleFun1 funDef) {
		
	}
	
	/**
	 * Construct a math function with one variable and it's derivative
	 * 
	 * @param funDef
	 * @param derivativeDef
	 */
	public MathFunTODO(SimpleFun1 funDef,
			SimpleFun1 derivativeDef) {
	}	
	
	
	public MathFunTODO(SimpleFun2 funDef) {
		
	}
	
	public MathFunTODO(SimpleFun2 funDef,
			SimpleFun2 derivativeX, 
			SimpleFun2 derivativeY) {
		
	}
	
	public MathFunTODO(SimpleFun3 funDef) {
		
	}
	
	public MathFunTODO(SimpleFun3 funDef,
			SimpleFun3 derivativeX,
			SimpleFun3 derivativeY,
			SimpleFun3 derivativeZ) {
		
	}

	@Override
	public double apply(Variable v) {
		// TODO Auto-generated method stub
		return 0;
	}
	
	public double value(double ...vars) {
		return 0;
	}
}
