package edu.uta.futureye.core;

import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.function.intf.ScalarShapeFunction;
import edu.uta.futureye.function.operator.FMath;
import edu.uta.futureye.lib.shapefun.SFBilinearLocal2D;
import edu.uta.futureye.lib.shapefun.SFLinearLocal1D;
import edu.uta.futureye.lib.shapefun.SFLinearLocal2D;
import edu.uta.futureye.lib.shapefun.SFLinearLocal3D;
import edu.uta.futureye.util.FutureyeException;
import edu.uta.futureye.util.list.VertexList;

public class CoordinateTransform {
	private List<String> fromVarNames = null;
	private List<String> toVarNames = null;
	private List<Function> transFuns = null;
	
	private ScalarShapeFunction sf1d1 = null;
	private ScalarShapeFunction sf1d2 = null;

	private ScalarShapeFunction sf2d1 = null;
	private ScalarShapeFunction sf2d2 = null;
	private ScalarShapeFunction sf2d3 = null;
	
	private ScalarShapeFunction sfb2d1 = null;
	private ScalarShapeFunction sfb2d2 = null;
	private ScalarShapeFunction sfb2d3 = null;
	private ScalarShapeFunction sfb2d4 = null;

	private ScalarShapeFunction sf3d1 = null;
	private ScalarShapeFunction sf3d2 = null;
	private ScalarShapeFunction sf3d3 = null;
	private ScalarShapeFunction sf3d4 = null;

	/**
	 * 构造一个坐标变换：
	 * e.g. 
	 *  2D (x,y)   -> (r,s)
	 *  3D (x,y,z) -> (r,s,t)
	 * @param fromVarName 
	 * @param toVarName
	 */
	public CoordinateTransform(List<String> fromVarNames, List<String> toVarNames) {
		if(fromVarNames.size() != toVarNames.size()) {
			FutureyeException e = new FutureyeException("ERROR: Number of variables mismatch.");
			e.printStackTrace();
			return;
		}
		this.fromVarNames = fromVarNames;
		this.toVarNames = toVarNames;
	}
	
	public CoordinateTransform(String[] fromVarNames, String[] toVarNames) {
		List<String> f = new LinkedList<String>();
		for(String s : fromVarNames)
			f.add(s);
		List<String> t = new LinkedList<String>();
		for(String s : toVarNames)
			t.add(s);
		
		this.fromVarNames = f;
		this.toVarNames = t;
	}	
	
	public CoordinateTransform(int dim) {
		fromVarNames = new LinkedList<String>();
		toVarNames = new LinkedList<String>();
		String[] f = {"x","y","z"};
		String[] t = {"r","s","t"};
		for(int i=0;i<dim;i++) {
			fromVarNames.add(f[i]);
			toVarNames.add(t[i]);
		}
	}

	/**
	 * 一维线性坐标变换
	 * @param be
	 * @return
	 */
	public Map<Vertex,ScalarShapeFunction> getTransformLinear1DShapeFunction(Element be) {
		Map<Vertex, ScalarShapeFunction> mapVS = new LinkedHashMap<Vertex, ScalarShapeFunction>();
		VertexList vl = be.vertices();
		if(vl.size() == 2) {
			if(sf1d1 == null) sf1d1 = new SFLinearLocal1D(1);
			if(sf1d2 == null) sf1d2 = new SFLinearLocal1D(2);
			sf1d1.asignElement(be);
			sf1d2.asignElement(be);
			mapVS.put(vl.at(1), sf1d1);
			mapVS.put(vl.at(2), sf1d2);
		}
		else {
			System.out.println("ERROR: transformLinear1D(), vl.size() != 2");
		}
		return mapVS;
	}
	
	/**
	 * 二维单元上线性坐标变换
	 * @param Element e
	 * @return 坐标变换所需的 “结点<=>形函数”对
	 */
	public Map<Vertex,ScalarShapeFunction> getTransformLinear2DShapeFunction(Element e) {
		VertexList vl = e.vertices();
		Map<Vertex, ScalarShapeFunction> mapVS = new LinkedHashMap<Vertex, ScalarShapeFunction>();
		if(vl.size() == 3) {
			if(sf2d1 == null) sf2d1 = new SFLinearLocal2D(1);
			if(sf2d2 == null) sf2d2 = new SFLinearLocal2D(2);
			if(sf2d3 == null) sf2d3 = new SFLinearLocal2D(3);
			sf2d1.asignElement(e);
			sf2d2.asignElement(e);
			sf2d3.asignElement(e);
			mapVS.put(vl.at(1), sf2d1);
			mapVS.put(vl.at(2), sf2d2);
			mapVS.put(vl.at(3), sf2d3);
		} else if(vl.size() == 4) {
			if(sfb2d1 == null) sfb2d1 = new SFBilinearLocal2D(1);
			if(sfb2d2 == null) sfb2d2 = new SFBilinearLocal2D(2);
			if(sfb2d3 == null) sfb2d3 = new SFBilinearLocal2D(3);
			if(sfb2d4 == null) sfb2d4 = new SFBilinearLocal2D(4);
			sfb2d1.asignElement(e);
			sfb2d2.asignElement(e);
			sfb2d3.asignElement(e);
			sfb2d4.asignElement(e);
			mapVS.put(vl.at(1), sfb2d1);
			mapVS.put(vl.at(2), sfb2d2);
			mapVS.put(vl.at(3), sfb2d3);
			mapVS.put(vl.at(4), sfb2d4);			
		} else if(vl.size() == 5) {
			//TODO
			System.out.println("ERROR: getTransformLinear2DShapeFunction() vl.size() == "+vl.size());
		} else {
			FutureyeException ex = new FutureyeException("ERROR: getTransformLinear2DShapeFunction() vl.size() == "+vl.size());
			ex.printStackTrace();
			System.exit(0);
		}
		return mapVS;
	}
	
	/**
	 * 三维单元上线性坐标变换
	 * @param Element e
	 * @return 坐标变换所需的 “结点<=>形函数”对
	 */
	public Map<Vertex,ScalarShapeFunction> getTransformLinear3DShapeFunction(Element e) {
		VertexList vl = e.vertices();
		Map<Vertex, ScalarShapeFunction> mapVS = new LinkedHashMap<Vertex, ScalarShapeFunction>();
		if(vl.size() == 4) {
			if(sf3d1 == null) sf3d1 = new SFLinearLocal3D(1);
			if(sf3d2 == null) sf3d2 = new SFLinearLocal3D(2);
			if(sf3d3 == null) sf3d3 = new SFLinearLocal3D(3);
			if(sf3d4 == null) sf3d4 = new SFLinearLocal3D(4);
			//TODO 可以省略？
			sf3d1.asignElement(e);
			sf3d2.asignElement(e);
			sf3d3.asignElement(e);
			sf3d4.asignElement(e);
			////////////////////////
			mapVS.put(vl.at(1), sf3d1);
			mapVS.put(vl.at(2), sf3d2);
			mapVS.put(vl.at(3), sf3d3);
			mapVS.put(vl.at(4), sf3d4);
		} else {
			FutureyeException ex = new FutureyeException("Error: getTransformLinear3DShapeFunction");
			ex.printStackTrace();
			System.exit(0);
		}
		return mapVS;
	}
	
	/**
	 * 利用单元上的形函数进行坐标变换
	 * @param e
	 * @return
	 */
	public Map<Vertex,ScalarShapeFunction> getTransformShapeFunctionByElement(Element e) {
		Map<Vertex, ScalarShapeFunction> mapVS = new LinkedHashMap<Vertex, ScalarShapeFunction>();
		VertexList vertices = e.vertices();
		for(int i=1;i<=vertices.size();i++) {
			//TODO 用顶点结点上的第一个自由度进行坐标变换
			mapVS.put(
					vertices.at(i), 
					e.getNodeDOFList(i).at(1).getSSF()
					);
		}
		return mapVS;
	}
	
	/**
	 * e.g. 二维笛卡尔坐标(x,y) -> 三角形面积（重心）坐标(r,s,t)
	 * (xi,yi), i=1,2,3 三角形三个顶点
	 * x = x(r,s,t) = x1*sf1 + x2*sf2 + x3*sf3 = x1*r + x2*s + x3*t
	 * y = y(r,s,t) = y1*sf1 + y2*sf2 + y3*sf3 = y1*r + y2*s + y3*t
	 * 利用后面的等式，可以简化运算
	 * 
	 * e.g. 三维笛卡尔坐标(x,y,z) -> 四面体体积（重心）坐标(r,s,t,u)
	 * (xi,yi), i=1,2,3,4 四面体四个顶点
	 * x = x(r,s,t,u) = x1*sf1 + x2*sf2 + x3*sf3 * x4*sf4 = x1*r + x2*s + x3*t * x4*u 
	 * y = y(r,s,t,u) = y1*sf1 + y2*sf2 + y3*sf3 * y4*sf4 = y1*r + y2*s + y3*t * y4*u 
	 * z = z(r,s,t,u) = z1*sf1 + z2*sf2 + z3*sf3 * z4*sf4 = z1*r + z2*s + z3*t * z4*u 
	 * 利用后面的等式，可以简化运算
	 * 
	 * @param mapVS
	 * @return 坐标变换函数列表
	 */
	public List<Function> getTransformFunction(Map<Vertex,ScalarShapeFunction> mapVS) {
		List<Function> r = new LinkedList<Function>();
		int dim = fromVarNames.size();
		for(int i=1;i<=dim;i++) {
			Function a = new FC(0.0);
			for(Entry<Vertex,ScalarShapeFunction> e : mapVS.entrySet()) {
				Vertex v = e.getKey();
				ScalarShapeFunction shapeFun = e.getValue();
				Function fm = FMath.Mult(
							new FC(v.coord(i)), shapeFun
						);
				a = FMath.Plus(a, fm);
			}
			//???? a.setVarNames(toVarNames);
			r.add(a);
		}
		return r;	
	}	
	
	public void setTransformFunction(List<Function> trans) {
		this.transFuns = trans;
	}
	
	public void transformLinear1D(Element be) {
		setTransformFunction(
				getTransformFunction(
					getTransformLinear1DShapeFunction(be)
					)
				);
	}
	
	public void transformLinear2D(Element e) {
		setTransformFunction(
				getTransformFunction(
					getTransformLinear2DShapeFunction(e)
					)
				);
	}
	public void transformLinear3D(Element e) {
		setTransformFunction(
				getTransformFunction(
					getTransformLinear3DShapeFunction(e)
					)
				);
	}
	
	public Function getJacobian1D() {
		Function[] funs = new Function[2];
		int index = 0;
		if(transFuns.size() != 2)
			return null;
		/*
		 * x = x(r,s,t) = x1*sf1 + x2*sf2 + x3*sf3 = x1*r + x2*s + x3*t
		 * y = y(r,s,t) = y1*sf1 + y2*sf2 + y3*sf3 = y1*r + y2*s + y3*t
		 */
		for(Function transFun : transFuns) {
			for(String var : toVarNames) { //free variable r,s
				funs[index++] = transFun.d(var);
			}
		}
		/*  
		 *  det(Jac)= Sqrt(0^2 + 1^2)
		 */
		return FMath.Sqrt(FMath.Plus(
				FMath.Mult(funs[0], funs[0]),
				FMath.Mult(funs[1], funs[1]))
				);
	}
	
	/**
	 * Dependent on setTransformFunction()
	 * @return
	 */
	public Function getJacobian2D() {
		Function[] funs = new Function[4];
		int index = 0;
		if(transFuns==null || transFuns.size() != 2)
			return null;
		for(Function transFun : transFuns) {
			for(String var : toVarNames) {
				funs[index++] = transFun.d(var);
			}
		}
		/*
		 * 要求结点编号为逆时针：funs[0:3]
		 *            |0 1|
		 *  det(Jac)= |2 3| = 0*3-1*2
		 */
		return FMath.Minus(
				FMath.Mult(funs[0], funs[3]),
				FMath.Mult(funs[1], funs[2]));
	}
	
	
	
	/**
	 * Dependent on setTransformFunction()
	 * @return
	 */
	public Function getJacobian3D() {
		Function[] funs = new Function[9];
		int index = 0;
		if(transFuns==null || transFuns.size() != 3)
			return null;
		/*
		 * x = x(r,s,t,u) = x1*sf1 + x2*sf2 + x3*sf3 * x4*sf4 = x1*r + x2*s + x3*t * x4*u 
		 * y = y(r,s,t,u) = y1*sf1 + y2*sf2 + y3*sf3 * y4*sf4 = y1*r + y2*s + y3*t * y4*u 
		 * z = z(r,s,t,u) = z1*sf1 + z2*sf2 + z3*sf3 * z4*sf4 = z1*r + z2*s + z3*t * z4*u 
		 */
		 for(Function transFun : transFuns) {
			for(String var : toVarNames) { //free variable r,s,t
				funs[index++] = transFun.d(var);
			}
		}
		/*
		 * 要求结点编号为逆时针 ：funs[0:8]
		 *            |0 1 2|
		 *  det(Jac)= |3 4 5| = 0*3-1*2
		 *            |6 7 8|
		 */
		Function det4875 = FMath.Minus(
				FMath.Mult(funs[4], funs[8]),
				FMath.Mult(funs[7], funs[5]));
		Function det6538 = FMath.Minus(
				FMath.Mult(funs[6], funs[5]),
				FMath.Mult(funs[3], funs[8]));
		Function det3764 = FMath.Minus(
				FMath.Mult(funs[3], funs[7]),
				FMath.Mult(funs[6], funs[4]));

		return FMath.Plus(FMath.Mult(funs[0], det4875), 
			   FMath.Plus(FMath.Mult(funs[1], det6538), 
					                 FMath.Mult(funs[2], det3764)
					            ));
	}
	
	public Function getJacobian3DFast(Element e) {
		double x1,x2,x3,x4;
		double y1,y2,y3,y4;
		double z1,z2,z3,z4;
		x1 = e.nodes.at(1).coord(1);
		x2 = e.nodes.at(2).coord(1);
		x3 = e.nodes.at(3).coord(1);
		x4 = e.nodes.at(4).coord(1);
		y1 = e.nodes.at(1).coord(2);
		y2 = e.nodes.at(2).coord(2);
		y3 = e.nodes.at(3).coord(2);
		y4 = e.nodes.at(4).coord(2);
		z1 = e.nodes.at(1).coord(3);
		z2 = e.nodes.at(2).coord(3);
		z3 = e.nodes.at(3).coord(3);
		z4 = e.nodes.at(4).coord(3);
		double[] items = new double[9];
		/*
		 * x = x(r,s,t,u) =  x1*r + x2*s + x3*t * x4*u 
		 * y = y(r,s,t,u) =  y1*r + y2*s + y3*t * y4*u 
		 * z = z(r,s,t,u) =  z1*r + z2*s + z3*t * z4*u 
		 */
		items[0] = x1 - x4;
		items[1] = x2 - x4;
		items[2] = x3 - x4;
		items[3] = y1 - y4;
		items[4] = y2 - y4;
		items[5] = y3 - y4;
		items[6] = z1 - z4;
		items[7] = z2 - z4;
		items[8] = z3 - z4;
		
		/*
		 * 要求结点编号为逆时针 ：items[0:8]
		 *            |0 1 2|
		 *  det(Jac)= |3 4 5| = 0*3-1*2
		 *            |6 7 8|
		 */
		double det4875 = items[4]*items[8]-items[7]*items[5];
		double det6538 = items[6]*items[5]-items[3]*items[8];
		double det3764 = items[3]*items[7]-items[6]*items[4];

		return new FC(
				  items[0]*det4875
				+ items[1]*det6538
				+ items[2]*det3764);
	}
}
