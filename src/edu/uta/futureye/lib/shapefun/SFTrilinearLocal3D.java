package edu.uta.futureye.lib.shapefun;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import edu.uta.futureye.core.CoordinateTransform;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.function.MultiVarFunc;
import edu.uta.futureye.function.VN;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.VariableArray;
import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.function.intf.ScalarShapeFunction;
import edu.uta.futureye.function.intf.ShapeFunction;
import edu.uta.futureye.util.FutureyeException;
import edu.uta.futureye.util.Utils;
import edu.uta.futureye.util.container.ObjList;

/**
 * Shape function for Hexahedron Element
 * 三线性形函数，用于六面体单元
 * 
 * 
 * @author liuyueming
 *
 */
public class SFTrilinearLocal3D extends MultiVarFunc implements ScalarShapeFunction {
	private int funIndex;
	private MathFunc funCompose = null;
	private MathFunc funOuter = null;
	private ObjList<String> innerVarNames = null;
	private double coef = 1.0;

	protected CoordinateTransform trans = new CoordinateTransform(3);
	public final static double[][] vt = {
			{ 1, 1, 1},
			{ 1,-1, 1},
			{-1,-1, 1},
			{-1, 1, 1},
			{ 1, 1,-1},
			{ 1,-1,-1},
			{-1,-1,-1},
			{-1, 1,-1}
		};
	MathFunc x_r = null;
	MathFunc x_s = null;
	MathFunc x_t = null;
	MathFunc y_r = null;
	MathFunc y_s = null;
	MathFunc y_t = null;
	MathFunc z_r = null;
	MathFunc z_s = null;
	MathFunc z_t = null;
	
	class InvJ extends MultiVarFunc {
		String rst,xyz;
		InvJ(String rst, String xyz) {
			this.rst = rst;
			this.xyz = xyz;
		}
		@Override
		public double apply(Variable v) {
			return apply(v,null);
		}
		@Override
		public double apply(Variable v, Map<Object,Object> cache) {
			Double detJ = null;
			double[][] J= null;
			if(cache != null) {
				detJ = (Double)cache.get(1);
				J = (double[][])cache.get(2);
			}
			if(detJ == null || J == null) {
				J = new double[3][3];
				J[0][0] = x_r.apply(v);
				J[0][1] = x_s.apply(v);
				J[0][2] = x_t.apply(v);
				J[1][0] = y_r.apply(v);
				J[1][1] = y_s.apply(v);
				J[1][2] = y_t.apply(v);
				J[2][0] = z_r.apply(v);
				J[2][1] = z_s.apply(v);
				J[2][2] = z_t.apply(v);	
				//@see ./doc/invA33.png
				detJ = Utils.determinant(J);
				if(cache != null) {
					cache.put(1, detJ);
					cache.put(2, J);
				}
			}

			if("r".equals(rst)) {
				if("x".equals(xyz))
					return (J[1][1]*J[2][2]-J[1][2]*J[2][1])/detJ;
				else if("y".equals(xyz))
					return (J[0][2]*J[2][1]-J[0][1]*J[2][2])/detJ;
				else if("z".equals(xyz))
					return (J[0][1]*J[1][2]-J[0][2]*J[1][1])/detJ;
			} else if("s".equals(rst)) {
				if("x".equals(xyz))
					return (J[1][2]*J[2][0]-J[1][0]*J[2][2])/detJ;
				else if("y".equals(xyz))
					return (J[0][0]*J[2][2]-J[0][2]*J[2][0])/detJ;
				else if("z".equals(xyz))
					return (J[0][2]*J[1][0]-J[0][0]*J[1][2])/detJ;
			} else {
				if("x".equals(xyz))
					return (J[1][0]*J[2][1]-J[1][1]*J[2][0])/detJ;
				else if("y".equals(xyz))
					return (J[0][1]*J[2][0]-J[0][0]*J[2][1])/detJ;
				else if("z".equals(xyz))
					return (J[0][0]*J[1][1]-J[0][1]*J[1][0])/detJ;
			}
			throw new FutureyeException("");
		}
		
		@Override
		public double[] applyAll(VariableArray valAry, Map<Object,Object> cache) {
			double[] detJ = null;
			double[][][] J= null;
			int len = valAry.length();
			double[] rlt = new double[len];
			if(cache != null) {
				detJ = (double[])cache.get(1);
				J = (double[][][])cache.get(2);
			}
			if(detJ == null || J == null) {
				J = new double[3][3][len];
				J[0][0] = x_r.applyAll(valAry,cache);
				J[0][1] = x_s.applyAll(valAry,cache);
				J[0][2] = x_t.applyAll(valAry,cache);
				J[1][0] = y_r.applyAll(valAry,cache);
				J[1][1] = y_s.applyAll(valAry,cache);
				J[1][2] = y_t.applyAll(valAry,cache);
				J[2][0] = z_r.applyAll(valAry,cache);
				J[2][1] = z_s.applyAll(valAry,cache);
				J[2][2] = z_t.applyAll(valAry,cache);
				//@see ./doc/invA33.png
				detJ = Utils.determinant(J);
				if(cache != null) {
					cache.put(1, detJ);
					cache.put(2, J);
				}
			}

			if("r".equals(rst)) {
				if("x".equals(xyz))
					for(int i=0;i<len;i++) rlt[i] = (J[1][1][i]*J[2][2][i]-J[1][2][i]*J[2][1][i])/detJ[i];
				else if("y".equals(xyz))
					for(int i=0;i<len;i++) rlt[i] = (J[0][2][i]*J[2][1][i]-J[0][1][i]*J[2][2][i])/detJ[i];
				else if("z".equals(xyz))
					for(int i=0;i<len;i++) rlt[i] = (J[0][1][i]*J[1][2][i]-J[0][2][i]*J[1][1][i])/detJ[i];
			} else if("s".equals(rst)) {
				if("x".equals(xyz))
					for(int i=0;i<len;i++) rlt[i] = (J[1][2][i]*J[2][0][i]-J[1][0][i]*J[2][2][i])/detJ[i];
				else if("y".equals(xyz))
					for(int i=0;i<len;i++) rlt[i] = (J[0][0][i]*J[2][2][i]-J[0][2][i]*J[2][0][i])/detJ[i];
				else if("z".equals(xyz))
					for(int i=0;i<len;i++) rlt[i] = (J[0][2][i]*J[1][0][i]-J[0][0][i]*J[1][2][i])/detJ[i];
			} else if("t".equals(rst)) {
				if("x".equals(xyz))
					for(int i=0;i<len;i++) rlt[i] = (J[1][0][i]*J[2][1][i]-J[1][1][i]*J[2][0][i])/detJ[i];
				else if("y".equals(xyz))
					for(int i=0;i<len;i++) rlt[i] = (J[0][1][i]*J[2][0][i]-J[0][0][i]*J[2][1][i])/detJ[i];
				else if("z".equals(xyz))
					for(int i=0;i<len;i++) rlt[i] = (J[0][0][i]*J[1][1][i]-J[0][1][i]*J[1][0][i])/detJ[i];
			} else {
				throw new FutureyeException();
			}
			return rlt;
		}
		@Override
		public int getOpOrder() {
			return OP_ORDER1;
		}
		@Override
		public String toString() {
			return rst+"_"+xyz;
		}
	}
	
	static class FOuter extends MultiVarFunc {
		int funIndex;
		FOuter(List<String> varNames, int funIndex) {
			this.varNames = varNames;
			this.funIndex = funIndex;
		}
		
		@Override
		public double apply(Variable v) {
			double r = v.get(VN.r);
			double s = v.get(VN.s);
			double t = v.get(VN.t);
			
			return (vt[funIndex][0]*r+1.0)*(vt[funIndex][1]*s+1.0)*(vt[funIndex][2]*t+1.0)/8.0;
		}
		
		@Override
		public double[] applyAll(VariableArray valAry, Map<Object,Object> cache) {
			int len = valAry.length();
			double[] r = valAry.get("r");
			double[] s = valAry.get("s");
			double[] t = valAry.get("t");
			double[] rlt = new double[len];
			double cr = vt[funIndex][0];
			double cs = vt[funIndex][1];
			double ct = vt[funIndex][2];
			for(int i=0;i<len;i++)
				rlt[i] = (cr*r[i]+1.0)*(cs*s[i]+1.0)*(ct*t[i]+1.0)/8.0;
			return rlt;
		}
		
		@Override
		public MathFunc diff(String var) {
			return new FOuter_d(varNames,var,funIndex);
		}
		public String toString() {
//			double cr = vt[funIndex][0];
//			double cs = vt[funIndex][1];
//			double ct = vt[funIndex][2];			
//			return String.format("(%.1f*r+1.0)*(%.1f*s+1.0)*(%.1f*t+1.0)/8.0", 
//					cr,cs,ct);
			return String.format("N%d", funIndex+1);
		}
	}
	
	static class FOuter_d extends MultiVarFunc {
		int funIndex;
		String var;
		FOuter_d(List<String> varNames, String var, int funIndex) {
			this.funIndex = funIndex;
			this.var = var;
		}
		@Override
		public double apply(Variable v) {
			double r = v.get(VN.r);
			double s = v.get(VN.s);
			double t = v.get(VN.t);
			
			double cr = vt[funIndex][0];
			double cs = vt[funIndex][1];
			double ct = vt[funIndex][2];
			
			if("r".equals(var))      return cr*(cs*s+1.0)*(ct*t+1.0)/8.0;
			else if("s".equals(var)) return (cr*r+1.0)*cs*(ct*t+1.0)/8.0;
			else if("t".equals(var)) return (cr*r+1.0)*(cs*s+1.0)*ct/8.0;
			else throw new FutureyeException("");
		}
		
		@Override
		public double[] applyAll(VariableArray valAry, Map<Object,Object> cache) {
			int len = valAry.length();
			double[] r = valAry.get("r");
			double[] s = valAry.get("s");
			double[] t = valAry.get("t");
			double[] rlt = new double[len];
			double cr = vt[funIndex][0];
			double cs = vt[funIndex][1];
			double ct = vt[funIndex][2];
			if("r".equals(var))      for(int i=0;i<len;i++) rlt[i] = cr*(cs*s[i]+1.0)*(ct*t[i]+1.0)/8.0;
			else if("s".equals(var)) for(int i=0;i<len;i++) rlt[i] = (cr*r[i]+1.0)*cs*(ct*t[i]+1.0)/8.0;
			else if("t".equals(var)) for(int i=0;i<len;i++) rlt[i] = (cr*r[i]+1.0)*(cs*s[i]+1.0)*ct/8.0;
			else throw new FutureyeException();
			return rlt;
		}
		@Override
		public int getOpOrder() {
			return OP_ORDER1;
		}
		public String toString() {
//			double cr = vt[funIndex][0];
//			double cs = vt[funIndex][1];
//			double ct = vt[funIndex][2];
			if("r".equals(var))
//				return String.format("%.1f*(%.1f*s+1.0)*(%.1f*t+1.0)/8.0", 
//						cr,cs,ct);
				return String.format("N%d_r", funIndex+1);
			else if("s".equals(var))
//				return String.format("(%.1f*r+1.0)*%.1f*(%.1f*t+1.0)/8.0", 
//						cr,cs,ct);
				return String.format("N%d_s", funIndex+1);
			else if("t".equals(var))
//				return String.format("(%.1f*r+1.0)*(%.1f*s+1.0)*%.1f/8.0", 
//						cr,cs,ct);
				return String.format("N%d_t", funIndex+1);
			else
				throw new FutureyeException("");
		}
	}
	
	/**
	 * 构造下列形函数中的一个：
	 *   Ni = (1+r*ri)*(1+s*si)*(1+t*ti)/8
	 * where
	 *   (ri,si,ti),i=1,...,8 are vertices coordinate of the hexahedron
	 * 
	 * @param funID = 1,...,8
	 * 
	 */	
	public void Create(int funID,double coef) {
		funIndex = funID - 1;
		if(funID<1 || funID>8) {
			System.out.println("ERROR: funID should be 1,...,8.");
			return;
		}
		
		varNames.add("r");
		varNames.add("s");
		varNames.add("t");
		innerVarNames = new ObjList<String>("x","y","z");
		
		//复合函数
		Map<String, MathFunc> fInners = new HashMap<String, MathFunc>(4);
		
		for(final String varName : varNames) {
			fInners.put(varName, new MultiVarFunc(innerVarNames.toList()) {
				
				//r_x,r_y,r_z, s_x,s_y,s_z, t_x,t_y,t_z
				public MathFunc diff(String var) {
/**
f(x,y,z) = g(r,s,t)
f_x = g_r*r_x + g_s*s_x + g_t*t_x ---(1)
f_y = g_r*r_y + g_s*s_y + g_t*t_y ---(2)
f_z = g_r*r_z + g_s*s_z + g_t*t_z ---(3)

for (1) let f=x,f=y,f=z we get 3 equations, solve them:
(x_r x_s x_t)   (r_x)   (1)
(y_r y_s y_z) * (s_x) = (0)
(z_r z_s z_t)   (t_x)   (0)

similarly, for (2):
(x_r x_s x_t)   (r_y)   (0)
(y_r y_s y_z) * (s_y) = (1)
(z_r z_s z_t)   (t_y)   (0)

for (3):
(x_r x_s x_t)   (r_z)   (0)
(y_r y_s y_z) * (s_z) = (0)
(z_r z_s z_t)   (t_z)   (1)

        (x_r x_s x_t)
Let J = (y_r y_s y_z)
        (z_r z_s z_t)

from the above 9 equations, we have:
 (r_x r_y r_z)
 (s_x s_y s_z) = inv(J)
 (t_x t_y t_z)

*/
					return new InvJ(varName,var);
				}

				@Override
				public double apply(Variable v) {
					// TODO Auto-generated method stub
					return 0;
				}
			});
		}
		
//		funOuter = new FAxpb("r",vt[funIndex][0]/2.0,0.5).M(
//				   new FAxpb("s",vt[funIndex][1]/2.0,0.5)).M(
//				   new FAxpb("t",vt[funIndex][2]/2.0,0.5));
		//速度提高1倍
		funOuter = new FOuter(varNames,funIndex);

		//使用复合函数构造形函数
		this.coef = coef;
		funCompose = FC.c(this.coef).M(
					funOuter.compose(fInners)
				);
	}

	public SFTrilinearLocal3D(int funID,double coef) {
		Create(funID,coef);
	}
	
	public SFTrilinearLocal3D(int funID) {
		Create(funID,1.0);
	}

	public MathFunc diff(String varName) {
		return funCompose.diff(varName);
	}

	public double apply(Variable v) {
		return funCompose.apply(v);
	}
	
	@Override
	public double[] applyAll(VariableArray v, Map<Object,Object> cache) {
		return funCompose.applyAll(v,cache);
	}

	@Override
	public void assignElement(Element e) {
//		//Coordinate transform and Jacbian on element e
//		List<Function> funs = trans.getTransformFunction(
//				trans.getTransformShapeFunctionByElement(e)
//				);
//		trans.setTransformFunction(funs);
//		
//		Function fx = funs.get(0); //x=x(r,s,t)
//		Function fy = funs.get(1); //y=y(r,s,t)
//		Function fz = funs.get(2); //z=z(r,s,t)
//		
//		x_r = fx._d("r");
//		x_s = fx._d("s");
//		x_t = fx._d("t");
//		y_r = fy._d("r");
//		y_s = fy._d("s");
//		y_t = fy._d("t");
//		z_r = fz._d("r");
//		z_s = fz._d("s");
//		z_t = fz._d("t");
		
		//faster
		MathFunc []funs = e.getCoordTrans().getJacobianMatrix();
		
		x_r = funs[0];
		x_s = funs[1];
		x_t = funs[2];
		y_r = funs[3];
		y_s = funs[4];
		y_t = funs[5];
		z_r = funs[6];
		z_s = funs[7];
		z_t = funs[8];
	}

	@Override
	public int getOpOrder() {
		return OP_ORDER1;
	}
	
	public String toString() {
		if(this.coef < 1.0)
			return this.coef+"*"+funOuter.toString();
		else
			return funOuter.toString();
	}

	SFBilinearLocal2D[] faceSF = {
			new SFBilinearLocal2D(1),
			new SFBilinearLocal2D(2),
			new SFBilinearLocal2D(3),
			new SFBilinearLocal2D(4)
		};
	
	@Override
	public ShapeFunction restrictTo(int funID) {
		return faceSF[funID-1];
	}

	@Override
	public ObjList<String> innerVarNames() {
		return innerVarNames;
	}
}
