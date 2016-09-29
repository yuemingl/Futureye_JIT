package edu.uta.futureye.lib.shapefun;

import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.function.AbstractVectorFunc;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.SpaceVectorFunction;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.function.intf.ScalarShapeFunction;
import edu.uta.futureye.function.intf.ShapeFunction;
import edu.uta.futureye.function.intf.VectorFunction;
import edu.uta.futureye.function.intf.VectorShapeFunction;
import edu.uta.futureye.util.FutureyeException;
import edu.uta.futureye.util.container.ObjList;

/**
 * P2/P0 Element
 * -Continuous quadratic velocity
 * -Piecewise constant pressure
 * 
 * * Velocity: Quadratic shape function: SFQuadraticLocal2DFast
 * * 速度：三角形局部坐标，二次函数
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
 * * Pressure: Piecewise constant shape function: SFConstant1
 * * 压强：分片常数型函数
 * NP=1
 * 
 * * 2D vector valued shape functions
 * * 二维单元上的形函数，速度压强共13个自由度：
 * Ni = (v1,v2,p)', i=1,...,13
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
 * N13 =  (0, 0, NP)'
 *
 * @author liuyueming
 */
public class QuadraticV_ConstantP extends AbstractVectorFunc 
								implements VectorShapeFunction {
	//(u1,u2,p)
	protected SpaceVectorFunction sf = null;
	protected int funIndex;
	protected ObjList<String> innerVarNames = null;
	
	/**
	 * 
	 * @param funID
	 */
	public QuadraticV_ConstantP(int funID) {
		funIndex = funID - 1;
		dim = 3;
		sf = new SpaceVectorFunction(dim);
		
		varNames.add("r");
		varNames.add("s");
		varNames.add("t");
		innerVarNames = new ObjList<String>("x","y");

		if(funIndex>=0 && funIndex<=5) {
			sf.set(1, new SFQuadraticLocal2DFast(funIndex+1));
			sf.set(2, new SFConstant0(innerVarNames));
			sf.set(3, new SFConstant0(innerVarNames));
		} else if(funIndex>=6 && funIndex<=11) {
			sf.set(1, new SFConstant0(innerVarNames));
			sf.set(2, new SFQuadraticLocal2DFast(funIndex-5));
			sf.set(3, new SFConstant0(innerVarNames));
		} else if(funIndex==12) {
			sf.set(1, new SFConstant0(innerVarNames));
			sf.set(2, new SFConstant0(innerVarNames));
			sf.set(3, new SFConstant1(innerVarNames));
		} else {
			throw new FutureyeException("ERROR: funID should be 1,2,3...,13.");
		}
	}
	
	@Override
	public void assignElement(Element e) {
		if(funIndex>=0 && funIndex<=5) {
			((ScalarShapeFunction)sf.get(1)).assignElement(e);
		} else if(funIndex>=6 && funIndex<=11) {
			((ScalarShapeFunction)sf.get(2)).assignElement(e);
		} else if(funIndex==12) {
			((ScalarShapeFunction)sf.get(3)).assignElement(e);
		}
	}

	@Override
	public ObjList<String> innerVarNames() {
		return innerVarNames;
	}

	QuadraticV_ConstantP1D []sf1D = null;
	/**
	 * Restrict to boundary edge:
	 * * Velocity:
	 * | \
	 * |  \
	 * |   \
	 * |    \
	 * |     \
	 * 1--3-- 2
	 * NV1,NV2,NV3
	 * 
	 * * Pressure:
	 * NP=1
	 * 
	 * Ni = (v1,v2,p)', i=1,...,7, on boundary
	 * 
	 * N1 =  (NV1,0,0)'
	 * N2 =  (NV2,0,0)'
	 * N3 =  (NV3,0,0)'
	 * N4 =  (0,NV1,0)'
	 * N5 =  (0,NV2,0)'
	 * N6 =  (0,NV3,0)'
	 * N7 =  (0,0,NP)'
	 * 
	 * @param funIndex
	 * @return
	 */
	@Override
	public ShapeFunction restrictTo(int funIndex) {
		if(sf1D == null) {
			sf1D = new QuadraticV_ConstantP1D[7];
			for(int i=1;i<=7;i++)
				sf1D[i-1] = new QuadraticV_ConstantP1D(i);
		}
		return sf1D[funIndex-1];
	}
	
	/**
	 * 需要定义边界单元，因为二维情况，向量值形函数dim=3，如果直接使用一维的，dim=2
	 * 
	 * @author liuyueming
	 *
	 */
	public class QuadraticV_ConstantP1D extends AbstractVectorFunc 
						implements VectorShapeFunction {
		//(u1,u2,p)
		protected SpaceVectorFunction sf = new SpaceVectorFunction(3);
		protected int funIndex;
		protected ObjList<String> innerVarNames = null;
		
		public QuadraticV_ConstantP1D(int funID) {
			dim = 3;
			
			this.funIndex = funID - 1;
			varNames.add("r");
			innerVarNames = new ObjList<String>("x");
			
			if(funIndex>=0 && funIndex<=2) {
				sf.set(1, new SFQuadraticLocal1D(funIndex+1));
				sf.set(2, new SFConstant0(innerVarNames));
				sf.set(3, new SFConstant0(innerVarNames));
			} else if(funIndex>=3 && funIndex<=5) {
				sf.set(1, new SFConstant0(innerVarNames));
				sf.set(2, new SFQuadraticLocal1D(funIndex-2));
				sf.set(3, new SFConstant0(innerVarNames));
			} else if(funIndex>=6 && funIndex<=6) {
				sf.set(1, new SFConstant0(innerVarNames));
				sf.set(2, new SFConstant0(innerVarNames));
				sf.set(3, new SFConstant1(innerVarNames));
			} else {
				throw new FutureyeException("ERROR: funID should be 1,2,3...,7.");
			}
			
		}

		@Override
		public void assignElement(Element e) {
			if(funIndex>=0 && funIndex<=2) {
				((ScalarShapeFunction)sf.get(1)).assignElement(e);
			} else if(funIndex>=3 && funIndex<=5) {
				((ScalarShapeFunction)sf.get(2)).assignElement(e);
			} else if(funIndex==6) {
				((ScalarShapeFunction)sf.get(3)).assignElement(e);
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
	};

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
