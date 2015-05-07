package edu.uta.futureye.lib.shapefun;

import java.util.Map;

import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.function.AbstractVectorFunc;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.VariableArray;
import edu.uta.futureye.function.basic.SpaceVectorFunction;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.function.intf.ScalarShapeFunction;
import edu.uta.futureye.function.intf.ShapeFunction;
import edu.uta.futureye.function.intf.VectorFunction;
import edu.uta.futureye.function.intf.VectorShapeFunction;
import edu.uta.futureye.util.FutureyeException;
import edu.uta.futureye.util.container.ObjList;

/**
 * 3D T1/P0 Element
 * -Continuous trilinear velocity
 * -Piecewise constant pressure
 * 
 * * Velocity: Trilinear shape function: SFTrilinearLocal3D
 * * 速度：六面体局部坐标，三线性函数
 * 
 * NVi(r,s,t) = (1+r*ri)*(1+s*si)*(1+t*ti)/8
 * where
 * (ri,si,ti),i=1,...,8 are vertices coordinate of the hexahedron
 * 
 * * Pressure: Piecewise constant shape function: SFConstant1
 * * 压强：分片常数型函数
 * NP=1
 * 
 * * 3D vector valued shape functions
 * * 三维单元上的形函数，速度压强共25个自由度：
 * Ni = (u1,u2,u3,p)', i=1,...,25
 * N_k       =  (NVk,0  ,0  ,0 )', k=1,...,8
 * N_(8 +k)  =  (0  ,0  ,0  ,0 )', k=1,...,8
 * N_(16+k)  =  (0  ,0  ,NVk,0 )', k=1,...,8
 * N_25      =  (0  ,0  ,0  ,NP)'
 *
 * @author liuyueming
 */
public class TrilinearV_ConstantP extends AbstractVectorFunc 
								implements VectorShapeFunction {
	//(u1,u2,u3,p)
	protected SpaceVectorFunction sf = null;
	protected int funIndex;
	protected ObjList<String> innerVarNames = null;
	
	/**
	 * 
	 * @param funID
	 */
	public TrilinearV_ConstantP(int funID) {
		funIndex = funID - 1;
		dim = 4;
		sf = new SpaceVectorFunction(dim);

		varNames.add("r");
		varNames.add("s");
		varNames.add("t");
		innerVarNames = new ObjList<String>("x","y","z");

		if(funIndex>=0 && funIndex<=7) {
			sf.set(1, new SFTrilinearLocal3D(funIndex+1));
			sf.set(2, new SFConstant0(innerVarNames));
			sf.set(3, new SFConstant0(innerVarNames));
			sf.set(4, new SFConstant0(innerVarNames));
		} else if(funIndex>=8 && funIndex<=15) {
			sf.set(1, new SFConstant0(innerVarNames));
			sf.set(2, new SFTrilinearLocal3D(funIndex-7));
			sf.set(3, new SFConstant0(innerVarNames));
			sf.set(4, new SFConstant0(innerVarNames));
		} else if(funIndex>=16 && funIndex<=23) {
			sf.set(1, new SFConstant0(innerVarNames));
			sf.set(2, new SFConstant0(innerVarNames));
			sf.set(3, new SFTrilinearLocal3D(funIndex-15));
			sf.set(4, new SFConstant0(innerVarNames));
		} else if(funIndex==24) {
			sf.set(1, new SFConstant0(innerVarNames));
			sf.set(2, new SFConstant0(innerVarNames));
			sf.set(3, new SFConstant0(innerVarNames));
			sf.set(4, new SFConstant1(innerVarNames));
		} else {
			throw new FutureyeException("ERROR: funID should be 1,2,3...,25.");
		}
	}
	
	@Override
	public void assignElement(Element e) {
		if(funIndex>=0 && funIndex<=7) {
			((ScalarShapeFunction)sf.get(1)).assignElement(e);
		} else if(funIndex>=8 && funIndex<=15) {
			((ScalarShapeFunction)sf.get(2)).assignElement(e);
		} else if(funIndex>=16 && funIndex<=23) {
			((ScalarShapeFunction)sf.get(3)).assignElement(e);
		} else if(funIndex==24) {
			((ScalarShapeFunction)sf.get(4)).assignElement(e);
		}
	}

	BilinearV_ConstantP []sf2D = null;
	
	@Override
	public ShapeFunction restrictTo(int funIndex) {
		if(sf2D == null) {
			sf2D = new BilinearV_ConstantP[13];
			for(int i=1;i<=13;i++)
				sf2D[i-1] = new BilinearV_ConstantP(i);
		}
		return sf2D[funIndex-1];
	}
	
	@Override
	public ObjList<String> innerVarNames() {
		return innerVarNames;
	}

	/**
	 * Restrict to boundary face of 3D element:
	 * 2D Q1/P0 Element
	 * -Continuous bilinear velocity
	 * -Piecewise constant pressure
	 * 
	 * * Velocity: Bilinear shape function: SFBilinearLocal2D
	 * * 速度：四边形局部坐标，双线性函数
	 * 
	 * 4----3
	 * |    |
	 * |    |
	 * 1----2
	 * 
	 * NV = NV(r,s) = NV( r(x,y), s(x,y) )
	 * NV1 = (1-r)*(1-s)/4
	 * NV2 = (1+r)*(1-s)/4
	 * NV3 = (1+r)*(1+s)/4
	 * NV4 = (1-r)*(1+s)/4
	 * 
	 * * Pressure: Piecewise constant shape function: SFConstant1
	 * * 压强：分片常数型函数
	 * NP=1
	 * 
	 * * 2D vector valued shape functions
	 * * 三维单元限制到二维面上的形函数，速度压强共13个自由度：
	 * Ni = (u1,u2,p)', i=1,...,13
	 * 
	 * N1  =  (NV1,   0,   0, 0)'
	 * N2  =  (NV2,   0,   0, 0)'
	 * N3  =  (NV3,   0,   0, 0)'
	 * N4  =  (NV4,   0,   0, 0)'
	 * N5  =  (  0, NV1,   0, 0)'
	 * N6  =  (  0, NV2,   0, 0)'
	 * N7  =  (  0, NV3,   0, 0)'
	 * N8  =  (  0, NV4,   0, 0)'
	 * N9  =  (  0,   0, NV1, 0)'
	 * N10 =  (  0,   0, NV2, 0)'
	 * N11 =  (  0,   0, NV3, 0)'
	 * N12 =  (  0,   0, NV4, 0)'
	 * N13 =  (  0,   0,   0,NP)'
	 *
	 */
	public class BilinearV_ConstantP extends AbstractVectorFunc 
									implements VectorShapeFunction {
		//(u1,u2,u3,p) 3D单元的边界单元不能直接使用2D的形函数，因为形函数需要向量dim=4，而2D的dim=3 (u1,u2, p)
		protected SpaceVectorFunction sf = null;
		protected int funIndex;
		protected ObjList<String> innerVarNames = null;
		
		/**
		 * 
		 * @param funID
		 */
		public BilinearV_ConstantP(int funID) {
			funIndex = funID - 1;
			dim = 4;
			sf = new SpaceVectorFunction(dim);
			
			varNames.add("r");
			varNames.add("s");
			innerVarNames = new ObjList<String>("x","y");

			if(funIndex>=0 && funIndex<=3) {
				sf.set(1, new SFBilinearLocal2D(funIndex+1));
				sf.set(2, new SFConstant0(innerVarNames));
				sf.set(3, new SFConstant0(innerVarNames));
				sf.set(4, new SFConstant0(innerVarNames));
			} else if(funIndex>=4 && funIndex<=7) {
				sf.set(1, new SFConstant0(innerVarNames));
				sf.set(2, new SFBilinearLocal2D(funIndex-3));
				sf.set(3, new SFConstant0(innerVarNames));
				sf.set(4, new SFConstant0(innerVarNames));
			} else if(funIndex>=8 && funIndex<=11) {
				sf.set(1, new SFConstant0(innerVarNames));
				sf.set(2, new SFConstant0(innerVarNames));
				sf.set(3, new SFBilinearLocal2D(funIndex-7));
				sf.set(4, new SFConstant0(innerVarNames));
			} else if(funIndex==12) {
				sf.set(1, new SFConstant0(innerVarNames));
				sf.set(2, new SFConstant0(innerVarNames));
				sf.set(3, new SFConstant0(innerVarNames));
				sf.set(3, new SFConstant1(innerVarNames));
			} else {
				throw new FutureyeException("ERROR: funID should be 1,2,3...,13.");
			}
		}
		
		@Override
		public void assignElement(Element e) {
			if(funIndex>=0 && funIndex<=3) {
				((ScalarShapeFunction)sf.get(1)).assignElement(e);
			} else if(funIndex>=4 && funIndex<=7) {
				((ScalarShapeFunction)sf.get(2)).assignElement(e);
			} else if(funIndex>=8 && funIndex<=11) {
				((ScalarShapeFunction)sf.get(3)).assignElement(e);
			} else if(funIndex==12) {
				((ScalarShapeFunction)sf.get(4)).assignElement(e);
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
	
	@Override
	public Vector[] valueArray(VariableArray v, Map<Object,Object> cache) {
		return (Vector[]) this.sf.valueArray(v,cache);
	}
	
	public String toString() {
		return sf.toString();
	}
}
