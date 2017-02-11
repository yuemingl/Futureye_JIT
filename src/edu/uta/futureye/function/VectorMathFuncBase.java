package edu.uta.futureye.function;

import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import edu.uta.futureye.algebra.SpaceVector;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.basic.SpaceVectorFunction;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.function.intf.VectorMathFunc;
import edu.uta.futureye.util.FutureyeException;
import edu.uta.futureye.util.Utils;

/**
 * Abstract vector function implementation.
 * <p>
 * Notice: Some member functions such as algebra operations will return 
 * an instance of SpaceVectorFunction which extends from AbstractVectorFunction
 * 
 * @author liuyueming
 *
 */
public abstract class VectorMathFuncBase implements VectorMathFunc {
	protected int dim = 0;
	protected List<String> varNames = new LinkedList<String>();
	protected String fName = null;
	
	public VectorMathFuncBase() {
	}
	
	public VectorMathFuncBase(int dim) {
		this.dim = dim;
	}
	
	public VectorMathFuncBase(int dim, List<String> varNames) {
		this.dim = dim;
		this.varNames = varNames;
	}
	
	public VectorMathFuncBase(int dim, String varName, String ...aryVarNames) {
		this.dim = dim;
		varNames.add(varName);
		for(String s : aryVarNames)
			varNames.add(s);
	}
	
	@Override
	public Vector value(Variable v) {
		Vector rlt = new SpaceVector(dim);
		for(int i=1; i<=dim; i++)
			rlt.set(i, this.get(i).apply(v));
		return rlt;
	}
	
	@Override
	public Vector[] valueArray(VariableArray valAry, Map<Object,Object> cache) {
		int len = valAry.length();
		Vector[] rlt = new SpaceVector[len];
		for(int j=0;j<len;j++) rlt[j] = new SpaceVector(dim);
		for(int i=1; i<=dim; i++) {
			double [] vals = this.get(i).applyAll(valAry, cache);
			for(int j=0;j<len;j++) {
				rlt[j].set(i, vals[j]);
			}
		}
		return rlt;
	}
	
	@Override
	public int getDim() {
		return this.dim;
	}
	
	@Override
	public void setDim(int dim) {
		this.dim = dim;
	}
	
	@Override
	public void setVarNames(List<String> varNames) {
		this.varNames = varNames;
	}

	@Override
	public List<String> varNames() {
		return this.varNames;
	}
	
	@Override
	public VectorMathFunc compose(Map<String,MathFunc> fInners) {
		for(int i=1; i<=dim; i++)
			this.set(i,this.get(i).compose(fInners));
		return this;
	}
	
	/////////////////////////////////////////////
	@Override
	public MathFunc get(int index) {
		throw new FutureyeException("Component of AbstractVectorFunction is not define!");
	}
	
	@Override
	public VectorMathFunc set(VectorMathFunc v) {
		for(int i=1; i<=dim; i++)
			this.set(i,v.get(i));
		this.setVarNames(v.varNames());
		return this;
	}
	
	@Override
	public VectorMathFunc set(double a, VectorMathFunc v) {
		for(int i=1; i<=dim; i++) {
			this.set(i,v.get(i).M(a));
		}
		this.setVarNames(v.varNames());
		return this;
	}
	
	@Override
	public VectorMathFunc inc(VectorMathFunc g) {
		for(int i=1; i<=dim; i++) {
			this.set(i,this.get(i).A(g.get(i)));
		}
		varNames = Utils.mergeList(varNames, g.varNames());
		return this;
	}
	
	@Override
	public VectorMathFunc inc(double a, VectorMathFunc g) {
		for(int i=1; i<=dim; i++) {
			this.set(i,this.get(i).A(g.get(i).M(a)));
		}
		varNames = Utils.mergeList(varNames, g.varNames());
		return this;
	}
	
	@Override
	public VectorMathFunc scale(double a) {
		for(int i=1; i<=dim; i++) {
			this.set(i,this.get(i).M(a));
		}
		return this;
	}
	
	@Override
	public VectorMathFunc ax(double a) {
		for(int i=1; i<=dim; i++) {
			this.set(i,this.get(i).M(a));
		}
		return this;
	}
	
	@Override
	public VectorMathFunc axpy(double a, VectorMathFunc g) {
		for(int i=1; i<=dim; i++) {
			this.set(i,this.get(i).M(a).A(g.get(i)));
		}
		varNames = Utils.mergeList(varNames, g.varNames());
		return this;
	}
	
	@Override
	public MathFunc dot(VectorMathFunc b) {
		if(dim != b.getDim()) {
			Exception e = new FutureyeException("Dims between two vector functions must be same!");
			e.printStackTrace();
			return null;
		}
		MathFunc rlt = FMath.C0;
		for(int i=1;i<=dim;i++) {
			rlt = rlt.A(this.get(i).M(b.get(i)));
		}
		return rlt;
	}

	@Override
	public MathFunc dot(Vector b) {
		MathFunc rlt = FMath.C0;
		for(int i=1;i<=dim;i++) {
			rlt = rlt.A(this.get(i).M(b.get(i)));
		}
		return rlt;
	}
	
	/////////////////////////////////////////////
	
	@Override
	public VectorMathFunc A(VectorMathFunc g) {
		SpaceVectorFunction svf = new SpaceVectorFunction(dim);
		for(int i=1; i<=dim; i++) {
			svf.set(i, this.get(i).A(g.get(i)));
		}
		return svf;
	}
	
	@Override
	public VectorMathFunc A(Vector g) {
		SpaceVectorFunction svf = new SpaceVectorFunction(dim);
		for(int i=1; i<=dim; i++) {
			svf.set(i, this.get(i).A(g.get(i)));
		}
		return svf;
	}
	
	@Override
	public VectorMathFunc S(VectorMathFunc g) {
		SpaceVectorFunction svf = new SpaceVectorFunction(dim);
		for(int i=1; i<=dim; i++) {
			svf.set(i, this.get(i).S(g.get(i)));
		}
		return svf;
	}
	
	@Override
	public VectorMathFunc S(Vector g) {
		SpaceVectorFunction svf = new SpaceVectorFunction(dim);
		for(int i=1; i<=dim; i++) {
			svf.set(i, this.get(i).S(g.get(i)));
		}
		return svf;
	}
	

	@Override
	public VectorMathFunc M(VectorMathFunc g) {
		SpaceVectorFunction svf = new SpaceVectorFunction(dim);
		for(int i=1; i<=dim; i++) {
			svf.set(i, this.get(i).M(g.get(i)));
		}
		return svf;
	}
	@Override
	public VectorMathFunc M(Vector g) {
		SpaceVectorFunction svf = new SpaceVectorFunction(dim);
		for(int i=1; i<=dim; i++) {
			svf.set(i, this.get(i).M(g.get(i)));
		}
		return svf;		
	}
	

	@Override
	public VectorMathFunc D(VectorMathFunc g) {
		SpaceVectorFunction svf = new SpaceVectorFunction(dim);
		for(int i=1; i<=dim; i++) {
			svf.set(i, this.get(i).D(g.get(i)));
		}
		return svf;		
	}
	@Override
	public VectorMathFunc D(Vector g) {
		SpaceVectorFunction svf = new SpaceVectorFunction(dim);
		for(int i=1; i<=dim; i++) {
			svf.set(i, this.get(i).D(g.get(i)));
		}
		return svf;		
	}
	
	///////////////////////////////////////////////////
	
	@Override
	public VectorMathFunc copy() {
		throw new UnsupportedOperationException();
	}

	@Override
	public String getExpression() {
		String s = varNames.toString();
		String rlt = "[";
		for(int i=0;i<dim;i++) {
			rlt += "F("+s.substring(1, s.length()-1)+")";
			if(i<dim-1)
				rlt += ", ";
		}
		return rlt+"]";
	}
	
	@Override
	public String getFName() {
		return this.fName;
	}
	
	@Override
	public VectorMathFunc setFName(String name) {
		this.fName = name;
		return this;
	}
	
	@Override
	public String toString() {
		if(getFName() == null) {
			return getExpression();
		} else 
			return getFName();
	}
	//////////////Operator overloading support through Java-OO//////////////////

	public VectorMathFunc valueOf(int v) {
		SpaceVectorFunction f = new SpaceVectorFunction(1);
		f.set(1, new FC(v));
		return f;
	}
	public VectorMathFunc valueOf(long v) {
		SpaceVectorFunction f = new SpaceVectorFunction(1);
		f.set(1, new FC(v));
		return f;
	}
	public VectorMathFunc valueOf(float v) {
		SpaceVectorFunction f = new SpaceVectorFunction(1);
		f.set(1, new FC(v));
		return f;
	}
	public VectorMathFunc valueOf(double v) {
		SpaceVectorFunction f = new SpaceVectorFunction(1);
		f.set(1, new FC(v));
		return f;
	}

	public VectorMathFunc add(VectorMathFunc other) {
		return this.A(other);
	}
	public VectorMathFunc add(int other) {
		return this.A(new SpaceVectorFunction(new SpaceVector(this.dim).setAll(other)));
	}
	public VectorMathFunc addRev(int other) {
		return this.A(new SpaceVectorFunction(new SpaceVector(this.dim).setAll(other)));
	}
	public VectorMathFunc add(long other) {
		return this.A(new SpaceVectorFunction(new SpaceVector(this.dim).setAll(other)));
	}
	public VectorMathFunc addRev(long other) {
		return this.A(new SpaceVectorFunction(new SpaceVector(this.dim).setAll(other)));
	}	
	public VectorMathFunc add(float other) {
		return this.A(new SpaceVectorFunction(new SpaceVector(this.dim).setAll(other)));
	}
	public VectorMathFunc addRev(float other) {
		return this.A(new SpaceVectorFunction(new SpaceVector(this.dim).setAll(other)));
	}	
	public VectorMathFunc add(double other) {
		return this.A(new SpaceVectorFunction(new SpaceVector(this.dim).setAll(other)));
	}
	public VectorMathFunc addRev(double other) {
		return this.A(new SpaceVectorFunction(new SpaceVector(this.dim).setAll(other)));
	}
	
	public VectorMathFunc subtract(VectorMathFunc other) {
		return this.S(other);
	}
	public VectorMathFunc subtract(int other) {
		return this.S(new SpaceVectorFunction(new SpaceVector(this.dim).setAll(other)));
	}
	public VectorMathFunc subtractRev(int other) {
		return new SpaceVectorFunction(new SpaceVector(this.dim).setAll(other)).S(this);
	}
	public VectorMathFunc subtract(long other) {
		return this.S(new SpaceVectorFunction(new SpaceVector(this.dim).setAll(other)));
	}
	public VectorMathFunc subtractRev(long other) {
		return new SpaceVectorFunction(new SpaceVector(this.dim).setAll(other)).S(this);
	}	
	public VectorMathFunc subtract(float other) {
		return this.S(new SpaceVectorFunction(new SpaceVector(this.dim).setAll(other)));
	}
	public VectorMathFunc subtractRev(float other) {
		return new SpaceVectorFunction(new SpaceVector(this.dim).setAll(other)).S(this);
	}
	public VectorMathFunc subtract(double other) {
		return this.S(new SpaceVectorFunction(new SpaceVector(this.dim).setAll(other)));
	}
	public VectorMathFunc subtractRev(double other) {
		return new SpaceVectorFunction(new SpaceVector(this.dim).setAll(other)).S(this);
	}
	
	public VectorMathFunc multiply(VectorMathFunc other) {
		return this.M(other);
	}
	public VectorMathFunc multiply(int other) {
		return this.M(new SpaceVectorFunction(new SpaceVector(this.dim).setAll(other)));
	}
	public VectorMathFunc multiplyRev(int other) {
		return this.M(new SpaceVectorFunction(new SpaceVector(this.dim).setAll(other)));
	}
	public VectorMathFunc multiply(long other) {
		return this.M(new SpaceVectorFunction(new SpaceVector(this.dim).setAll(other)));
	}
	public VectorMathFunc multiplyRev(long other) {
		return this.M(new SpaceVectorFunction(new SpaceVector(this.dim).setAll(other)));
	}
	public VectorMathFunc multiply(float other) {
		return this.M(new SpaceVectorFunction(new SpaceVector(this.dim).setAll(other)));
	}
	public VectorMathFunc multiplyRev(float other) {
		return this.M(new SpaceVectorFunction(new SpaceVector(this.dim).setAll(other)));
	}
	public VectorMathFunc multiply(double other) {
		return this.M(new SpaceVectorFunction(new SpaceVector(this.dim).setAll(other)));
	}
	public VectorMathFunc multiplyRev(double other) {
		return this.M(new SpaceVectorFunction(new SpaceVector(this.dim).setAll(other)));
	}
	
	public VectorMathFunc divide(VectorMathFunc other) {
		return this.D(other);
	}	
	public VectorMathFunc divide(int other) {
		return this.D(new SpaceVectorFunction(new SpaceVector(this.dim).setAll(other)));
	}
	public VectorMathFunc divideRev(int other) {
		return new SpaceVectorFunction(new SpaceVector(this.dim).setAll(other)).D(this);
	}
	public VectorMathFunc divide(long other) {
		return this.D(new SpaceVectorFunction(new SpaceVector(this.dim).setAll(other)));
	}
	public VectorMathFunc divideRev(long other) {
		return new SpaceVectorFunction(new SpaceVector(this.dim).setAll(other)).D(this);
	}
	public VectorMathFunc divide(float other) {
		return this.D(new SpaceVectorFunction(new SpaceVector(this.dim).setAll(other)));
	}
	public VectorMathFunc divideRev(float other) {
		return new SpaceVectorFunction(new SpaceVector(this.dim).setAll(other)).D(this);
	}
	public VectorMathFunc divide(double other) {
		return this.D(new SpaceVectorFunction(new SpaceVector(this.dim).setAll(other)));
	}
	public VectorMathFunc divideRev(double other) {
		return new SpaceVectorFunction(new SpaceVector(this.dim).setAll(other)).D(this);
	}
	
	public VectorMathFunc negate() {
		return new SpaceVectorFunction(new SpaceVector(this.dim).setAll(0)).S(this);
	};
	
	/////////////////////////////////////////////////////////////
		
}
