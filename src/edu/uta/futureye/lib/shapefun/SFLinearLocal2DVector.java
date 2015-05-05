package edu.uta.futureye.lib.shapefun;

import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.function.AbstractVectorFunction;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.SpaceVectorFunction;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.function.intf.ScalarShapeFunction;
import edu.uta.futureye.function.intf.ShapeFunction;
import edu.uta.futureye.function.intf.VectorFunction;
import edu.uta.futureye.function.intf.VectorShapeFunction;
import edu.uta.futureye.util.container.ObjList;

public class SFLinearLocal2DVector extends AbstractVectorFunction 
								implements VectorShapeFunction {
	protected SpaceVectorFunction sf = new SpaceVectorFunction(2);
	protected int funIndex;
	protected ObjList<String> innerVarNames = null;

	
	public SFLinearLocal2DVector(int funID) {
		this.dim = 2;
		this.funIndex = funID - 1;
		
		varNames.add("r");
		varNames.add("s");
		innerVarNames = new ObjList<String>("x","y");

		if(funIndex>=0 && funIndex<=2) {
			sf.set(1, new SFLinearLocal2D(funIndex+1));
			sf.set(2, new SFConstant0(innerVarNames));
		} else if(funIndex>=3 && funIndex<=5) {
			sf.set(1, new SFConstant0(innerVarNames));
			sf.set(2, new SFLinearLocal2D(funIndex-2));
		}
		
	}
	
	@Override
	public void assignElement(Element e) {
		if(funIndex>=0 && funIndex<=2) {
			((ScalarShapeFunction)sf.get(1)).assignElement(e);
		} else if(funIndex>=3 && funIndex<=5) {
			((ScalarShapeFunction)sf.get(2)).assignElement(e);
		}
	}

	@Override
	public ObjList<String> innerVarNames() {
		return innerVarNames;
	}

	public class SFLinearLocal1DVector extends AbstractVectorFunction 
			implements VectorShapeFunction {
		protected SpaceVectorFunction sf = new SpaceVectorFunction(2);
		protected int funIndex;
		protected ObjList<String> innerVarNames = null;
		
		public SFLinearLocal1DVector(int funID) {
			this.dim = 2;
			this.funIndex = funID - 1;
			
			varNames.add("r");
			innerVarNames = new ObjList<String>("x");

			if(funIndex>=0 && funIndex<=1) {
				sf.set(1, new SFLinearLocal1D(funIndex+1));
				sf.set(2, new SFConstant0(innerVarNames));
			} else if(funIndex>=2 && funIndex<=3) {
				sf.set(1, new SFConstant0(innerVarNames));
				sf.set(2, new SFLinearLocal1D(funIndex-1));
			}
		}

		@Override
		public void assignElement(Element e) {
			if(funIndex>=0 && funIndex<=1) {
				((ScalarShapeFunction)sf.get(1)).assignElement(e);
			} else if(funIndex>=2 && funIndex<=3) {
				((ScalarShapeFunction)sf.get(2)).assignElement(e);
			}
		}
		@Override
		public ObjList<String> innerVarNames() {
			return innerVarNames;
		}
		
		@Override
		public Vector value(Variable v) {
			return (Vector) this.sf.value(v);
		}

		@Override
		public ShapeFunction restrictTo(int funIndex) {
			throw new UnsupportedOperationException();
		}
		
		@Override
		public MathFunc get(int index) {
			return sf.get(index);
		}

		@Override
		public void set(int index, MathFunc value) {
			throw new UnsupportedOperationException();
		}
		public String toString() {
			return sf.toString();
		}
		@Override
		public MathFunc dot(Vector b) {
			return sf.dot(b);
		}
	}
	
	SFLinearLocal1DVector []sf1D = null;
	
	@Override
	public ShapeFunction restrictTo(int funIndex) {
		if(sf1D == null) {
			sf1D = new SFLinearLocal1DVector[4];
			for(int i=1;i<=4;i++)
				sf1D[i-1] = new SFLinearLocal1DVector(i);
		}
		return sf1D[funIndex-1];
	}

	@Override
	public MathFunc dot(VectorFunction b) {
		return sf.dot(b);
	}

	@Override
	public MathFunc dot(Vector b) {
		return sf.dot(b);
	}

	@Override
	public MathFunc get(int index) {
		return sf.get(index);
	}

	@Override
	public void set(int index, MathFunc value) {
		throw new UnsupportedOperationException();
	}
	
	@Override
	public Vector value(Variable v) {
		return (Vector) this.sf.value(v);
	}
	
	public String toString() {
		return sf.toString();
	}

}
