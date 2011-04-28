package edu.uta.futureye.lib.shapefun;

import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.function.AbstractVectorFunction;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.SpaceVectorFunction;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.function.intf.ScalarShapeFunction;
import edu.uta.futureye.function.intf.ShapeFunction;
import edu.uta.futureye.function.intf.VectorFunction;
import edu.uta.futureye.function.intf.VectorShapeFunction;
import edu.uta.futureye.util.container.ObjList;

public class SFBilinearLocal2DVector extends AbstractVectorFunction 
								implements VectorShapeFunction {
	protected SpaceVectorFunction sf = new SpaceVectorFunction(2);
	protected int funIndex;
	protected ObjList<String> innerVarNames = null;

	
	public SFBilinearLocal2DVector(int funID) {
		this.dim = 2;
		this.funIndex = funID - 1;
		
		varNames.add("r");
		varNames.add("s");
		innerVarNames = new ObjList<String>("x","y");

		if(funIndex>=0 && funIndex<=3) {
			sf.set(1, new SFBilinearLocal2D(funIndex+1));
			sf.set(2, new SF0(innerVarNames));
		} else if(funIndex>=4 && funIndex<=7) {
			sf.set(1, new SF0(innerVarNames));
			sf.set(2, new SFBilinearLocal2D(funIndex-3));
		}
		
	}
	
	@Override
	public void asignElement(Element e) {
		if(funIndex>=0 && funIndex<=3) {
			((ScalarShapeFunction)sf.get(1)).asignElement(e);
		} else if(funIndex>=4 && funIndex<=7) {
			((ScalarShapeFunction)sf.get(2)).asignElement(e);
		}
	}

	@Override
	public ObjList<String> innerVarNames() {
		return innerVarNames;
	}

	public class SFBilinearLocal1DVector extends AbstractVectorFunction 
			implements VectorShapeFunction {
		protected SpaceVectorFunction sf = new SpaceVectorFunction(2);
		protected int funIndex;
		protected ObjList<String> innerVarNames = null;
		
		public SFBilinearLocal1DVector(int funID) {
			this.dim = 2;
			this.funIndex = funID - 1;
			
			varNames.add("r");
			innerVarNames = new ObjList<String>("x");

			if(funIndex>=0 && funIndex<=1) {
				sf.set(1, new SFLinearLocal1D(funIndex+1));
				sf.set(2, new SF0(innerVarNames));
			} else if(funIndex>=2 && funIndex<=3) {
				sf.set(1, new SF0(innerVarNames));
				sf.set(2, new SFLinearLocal1D(funIndex-1));
			}
		}

		@Override
		public void asignElement(Element e) {
			if(funIndex>=0 && funIndex<=1) {
				((ScalarShapeFunction)sf.get(1)).asignElement(e);
			} else if(funIndex>=2 && funIndex<=3) {
				((ScalarShapeFunction)sf.get(2)).asignElement(e);
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
		public Function get(int index) {
			return sf.get(index);
		}

		@Override
		public void set(int index, Function value) {
			throw new UnsupportedOperationException();
		}
		public String toString() {
			return sf.toString();
		}
		@Override
		public Function dot(Vector b) {
			return sf.dot(b);
		}
	}
	
	SFBilinearLocal1DVector []sf1D = null;
	
	@Override
	public ShapeFunction restrictTo(int funIndex) {
		if(sf1D == null) {
			sf1D = new SFBilinearLocal1DVector[4];
			for(int i=1;i<=4;i++)
				sf1D[i-1] = new SFBilinearLocal1DVector(i);
		}
		return sf1D[funIndex-1];
	}

	@Override
	public Function dot(VectorFunction b) {
		return sf.dot(b);
	}

	@Override
	public Function dot(Vector b) {
		return sf.dot(b);
	}

	@Override
	public Function get(int index) {
		return sf.get(index);
	}

	@Override
	public void set(int index, Function value) {
		throw new UnsupportedOperationException();
	}
	
	@Override
	public void print() {
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
