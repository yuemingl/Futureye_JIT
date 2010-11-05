package edu.uta.futureye.function.intf;

import java.util.List;

import edu.uta.futureye.function.Variable;

public interface Function {
	public double value(Variable v);
	
	public void setVarNames(List<String> varNames);
	public List<String> varNames();
}
