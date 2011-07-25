package edu.uta.futureye.function.basic;

import java.util.List;

import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.application.Tools;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.function.AbstractFunction;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.util.FutureyeException;

/**
 * Vector to Function
 * Evaluate function values based on vector indices in Variable v
 * 
 * 2011/6/27
 * + Function _d(String varName)
 * 
 * @author liuyueming
 *
 */
public class Vector2Function extends AbstractFunction {
	Vector u = null;
	Mesh mesh = null;
	
	public Vector2Function(Vector u) {
		this.u = u;
	}
	
	public Vector2Function(Vector u, Mesh mesh,
			String varName, String ...aryVarNames) {
		this.u = u;
		this.mesh = mesh;
		varNames.add(varName);
		for(String s : aryVarNames)
			varNames.add(s);
	}
	
	public Vector2Function(Vector u, Mesh mesh,
			List<String> varNames) {
		this.u = u;
		this.mesh = mesh;
		this.varNames = varNames;
	}	
	
	@Override
	public double value(Variable v) {
		int index = v.getIndex();
		if(index <= 0) {
			throw new FutureyeException("Error: Vector2Function index="+index);
		} else {
			return u.get(index);//注：下标错位会造成结果出现随机混乱
		}
	}
	
	@Override
	public Function _d(String varName) {
		Vector vd = Tools.computeDerivative(mesh, u, varName);
		Function fd = new Vector2Function(vd,mesh,this.varNames);
		return fd;
	}

}
