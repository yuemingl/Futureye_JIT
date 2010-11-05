package edu.uta.futureye.test;

import java.util.List;

import edu.uta.futureye.algebra.Vector;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.intf.Function;

public class VectorBasedFunction implements Function {
	Vector u = null;
	
	public VectorBasedFunction(Vector u) {
		this.u = u;
	}
	
	@Override
	public void setVarNames(List<String> varNames) {
		// TODO Auto-generated method stub

	}

	@Override
	public double value(Variable v) {
		int index = v.getIndex();
		if(index <= 0) {
			Exception e = new Exception("Error: VectorBasedFunction index="+index);
			e.printStackTrace();
			return 0;
		} else {
			//下标错位会造成结果出现随机混乱
			return u.get(index);
		}
	}

	@Override
	public List<String> varNames() {
		// TODO Auto-generated method stub
		return null;
	}
}
