package edu.uta.futureye.function.basic;

import java.util.HashMap;
import java.util.Map;

import edu.uta.futureye.function.AbstractFunction;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.util.FutureyeException;
import edu.uta.futureye.util.PairDoubleInteger;
import edu.uta.futureye.util.container.ObjList;

public class DiscreteIndexFunction extends AbstractFunction {
	protected Map<Integer, Double> data = new HashMap<Integer, Double>();
	
	public DiscreteIndexFunction() {
	}
	
	public DiscreteIndexFunction(ObjList<PairDoubleInteger> list) {
		for(int i=1;i<=list.size();i++) {
			PairDoubleInteger pair = list.at(i);
			this.data.put(pair.i,pair.d);
		}
	}
	
	public void set(int index, double value) {
		this.data.put(index, value);
	}
	
	@Override
	public double value(Variable v) {
		if(v.getIndex()<=0) 
			throw new FutureyeException("v.getIndex()="+v.getIndex());
		return data.get(v.getIndex());
	}

}
