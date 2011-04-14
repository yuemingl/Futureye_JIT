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

/**
 * Velocity: Quadratic shape function
 * 速度：三角形局部坐标，二次函数
 * 
 * 3
 * | \
 * |  \
 * 6   5
 * |    \
 * |     \
 * 1--4-- 2
 * 
 * NV = NV(r,s,t) = NV( r(x,y), s(x,y), t(x,y) )
 * NV1 = (2*r-1)*r
 * NV2 = (2*s-1)*s
 * NV3 = (2*t-1)*t
 * NV4 = 4*r*s
 * NV5 = 4*s*t
 * NV6 = 4*r*t
 * 
 * Pressure: Linear shape function
 * 压强：三角形局部坐标，线性型函数
 * 3
 * | \
 * |  \
 * |   \
 * |    \
 * 1---- 2
 * 
 * NP = NP(r,s,t) = NP( r(x,y), s(x,y), t(x,y) )
 * NP1 = r
 * NP2 = s
 * NP3 = t
 *
 * 2D vector valued shape functions
 * 二维单元上的形函数，速度压强共15个自由度：
 * Ni = (v1,v2,p)', i=1,...,15
 * 
 * N1  =  (NV1, 0, 0)'
 * N2  =  (NV2, 0, 0)'
 * N3  =  (NV3, 0, 0)'
 * N4  =  (NV4, 0, 0)'
 * N5  =  (NV5, 0, 0)'
 * N6  =  (NV6, 0, 0)'
 * N7  =  (0, NV1, 0)'
 * N8  =  (0, NV2, 0)'
 * N9  =  (0, NV3, 0)'
 * N10 =  (0, NV4, 0)'
 * N11 =  (0, NV5, 0)'
 * N12 =  (0, NV6, 0)'
 * N13 =  (0, 0, NP1)'
 * N14 =  (0, 0, NP2)'
 * N15 =  (0, 0, NP3)'
 *
 *
 * @author liuyueming
 */
public class QuadraticV_LinearP extends AbstractVectorFunction 
								implements VectorShapeFunction {
	//(u1,u2,p)
	protected SpaceVectorFunction sf = new SpaceVectorFunction(3);
	protected int funIndex;
	protected ObjList<String> innerVarNames = null;
	
	/**
	 * 
	 * @param funID
	 */
	public QuadraticV_LinearP(int funID) {
		dim = 3;
		
		funIndex = funID - 1;

		varNames.add("r");
		varNames.add("s");
		varNames.add("t");
		innerVarNames = new ObjList<String>("x","y");

		if(funIndex>=0 && funIndex<=5) {
			sf.set(1, new SFQuadraticLocal2DFast(funIndex+1));
			sf.set(2, new SF0(innerVarNames));
			sf.set(3, new SF0(innerVarNames));
		} else if(funIndex>=6 && funIndex<=11) {
			sf.set(1, new SF0(innerVarNames));
			sf.set(2, new SFQuadraticLocal2DFast(funIndex-5));
			sf.set(3, new SF0(innerVarNames));
		} else if(funIndex>=12 && funIndex<=14) {
			sf.set(1, new SF0(innerVarNames));
			sf.set(2, new SF0(innerVarNames));
			sf.set(3, new SFLinearLocal2D(funIndex-11));
		}


	}
	
	@Override
	public void asignElement(Element e) {
		if(funIndex>=0 && funIndex<=5) {
			((ScalarShapeFunction)sf.get(1)).asignElement(e);
		} else if(funIndex>=6 && funIndex<=11) {
			((ScalarShapeFunction)sf.get(2)).asignElement(e);
		} else if(funIndex>=12 && funIndex<=14) {
			((ScalarShapeFunction)sf.get(3)).asignElement(e);
		}
	}

	@Override
	public ObjList<String> innerVarNames() {
		return innerVarNames;
	}

	QuadraticV_LinearP1D []sf1D = null;
	/**
	 * Restrict to boundary
	 * | \
	 * |  \
	 * |   \
	 * |    \
	 * |     \
	 * 1--3-- 2
	 * NV1,NV2,NV3
	 * 
	 * | \
	 * |  \
	 * |   \
	 * |    \
	 * 1---- 2 
	 * NP1,NP2
	 * 
	 * Ni = (v1,v2,p)', i=1,...,8, on boundary
	 * 
	 * N1 =  (NV1,0,0)'
	 * N2 =  (NV2,0,0)'
	 * N3 =  (NV3,0,0)'
	 * N4 =  (0,NV1,0)'
	 * N5 =  (0,NV2,0)'
	 * N6 =  (0,NV3,0)'
	 * N7 =  (0,0,NP1)'
	 * N8 =  (0,0,NP2)'
	 * 
	 * @param funIndex
	 * @return
	 */
	@Override
	public ShapeFunction restrictTo(int funIndex) {
		if(sf1D == null) {
			sf1D = new QuadraticV_LinearP1D[8];
			for(int i=1;i<=8;i++)
				sf1D[i-1] = new QuadraticV_LinearP1D(i);
		}
		return sf1D[funIndex-1];
	}
	
	public class QuadraticV_LinearP1D extends AbstractVectorFunction implements VectorShapeFunction {
		//(u1,u2,p)
		protected SpaceVectorFunction sf = new SpaceVectorFunction(3);
		protected int funIndex;
		protected ObjList<String> innerVarNames = null;
		
		public QuadraticV_LinearP1D(int funID) {
			dim = 3;
			
			this.funIndex = funID - 1;
			varNames.add("r");
			innerVarNames = new ObjList<String>("x");
			
			if(funIndex>=0 && funIndex<=2) {
				sf.set(1, new SFQuadraticLocal1D(funIndex+1));
				sf.set(2, new SF0(innerVarNames));
				sf.set(3, new SF0(innerVarNames));
			} else if(funIndex>=3 && funIndex<=5) {
				sf.set(1, new SF0(innerVarNames));
				sf.set(2, new SFQuadraticLocal1D(funIndex-2));
				sf.set(3, new SF0(innerVarNames));
			} else if(funIndex>=6 && funIndex<=7) {
				sf.set(1, new SF0(innerVarNames));
				sf.set(2, new SF0(innerVarNames));
				sf.set(3, new SFLinearLocal1D(funIndex-5));
			}
			
		}

		@Override
		public void asignElement(Element e) {
			if(funIndex>=0 && funIndex<=2) {
				((ScalarShapeFunction)sf.get(1)).asignElement(e);
			} else if(funIndex>=3 && funIndex<=5) {
				((ScalarShapeFunction)sf.get(2)).asignElement(e);
			} else if(funIndex>=6 && funIndex<=7) {
				((ScalarShapeFunction)sf.get(3)).asignElement(e);
			}
		}
		
		@Override
		public Vector value(Variable v) {
			return (Vector) this.sf.value(v);
		}
		
		@Override
		public ObjList<String> innerVarNames() {
			return innerVarNames;
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
	};

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
