package edu.uta.futureye.function.basic;

import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.function.AbstractFunction;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.util.FutureyeException;

/**
 * Vector to Function, evaluate base on vector index in Variable v
 * 
 * @author liuyueming
 *
 */
public class Vector2Function extends AbstractFunction {
	Vector u = null;
	
	public Vector2Function(Vector u) {
		this.u = u;
	}
	
	@Override
	public double value(Variable v) {
		int index = v.getIndex();
		if(index <= 0) {
			FutureyeException ex = new FutureyeException("Error: Vector2Function index="+index);
			ex.printStackTrace();
			System.exit(-1);
		} else {
			return u.get(index);//注：下标错位会造成结果出现随机混乱
		}
		return 0.0;
	}
}
