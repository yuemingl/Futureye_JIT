package edu.uta.futureye.core;

import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import edu.uta.futureye.function.DerivativeIndicator;
import edu.uta.futureye.function.basic.FConstant;
import edu.uta.futureye.function.intf.FunctionDerivable;
import edu.uta.futureye.function.intf.ShapeFunction;
import edu.uta.futureye.function.operator.FOBasicDerivable;
import edu.uta.futureye.function.shape.SFBilinearLocal2D;
import edu.uta.futureye.function.shape.SFLinearLocal1D;
import edu.uta.futureye.function.shape.SFLinearLocal2D;

public class CoordinateTransform {
	private List<String> fromVarNames = null;
	private List<String> toVarNames = null;
	private List<FunctionDerivable> transFuns = null;
	
	private ShapeFunction sf1d1 = null;
	private ShapeFunction sf1d2 = null;

	private ShapeFunction sf2d1 = null;
	private ShapeFunction sf2d2 = null;
	private ShapeFunction sf2d3 = null;
	
	private ShapeFunction sfb2d1 = null;
	private ShapeFunction sfb2d2 = null;
	private ShapeFunction sfb2d3 = null;
	private ShapeFunction sfb2d4 = null;
	
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
			System.out.println("ERROR: Number of variables mismatch.");
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
	public Map<Vertex,ShapeFunction> getTransformLinear1DShapeFunction(Element be) {
		Map<Vertex, ShapeFunction> mapVS = new LinkedHashMap<Vertex, ShapeFunction>();
		List<Vertex> vl = be.getVertexList();
		if(vl.size() == 2) {
			if(sf1d1 == null) sf1d1 = new SFLinearLocal1D(1);
			if(sf1d2 == null) sf1d2 = new SFLinearLocal1D(2);
			sf1d1.asignElement(be);
			sf1d2.asignElement(be);
			mapVS.put(vl.get(0), sf1d1);
			mapVS.put(vl.get(1), sf1d2);
		}
		else {
			System.out.println("ERROR: transformLinear1D(), vl.size() != 2");
		}
		return mapVS;
	}
	
	/**
	 * 二维线性坐标变换
	 * @param e
	 * @return
	 */
	public Map<Vertex,ShapeFunction> getTransformLinear2DShapeFunction(Element e) {
		List<Vertex> vl = e.getVertexList();
		Map<Vertex, ShapeFunction> mapVS = new LinkedHashMap<Vertex, ShapeFunction>();
		if(vl.size() == 3) {
			if(sf2d1 == null) sf2d1 = new SFLinearLocal2D(1);
			if(sf2d2 == null) sf2d2 = new SFLinearLocal2D(2);
			if(sf2d3 == null) sf2d3 = new SFLinearLocal2D(3);
			sf2d1.asignElement(e);
			sf2d2.asignElement(e);
			sf2d3.asignElement(e);
			mapVS.put(vl.get(0), sf2d1);
			mapVS.put(vl.get(1), sf2d2);
			mapVS.put(vl.get(2), sf2d3);
		} else if(vl.size() == 4) {
			if(sfb2d1 == null) sfb2d1 = new SFBilinearLocal2D(1);
			if(sfb2d2 == null) sfb2d2 = new SFBilinearLocal2D(2);
			if(sfb2d3 == null) sfb2d3 = new SFBilinearLocal2D(3);
			if(sfb2d4 == null) sfb2d4 = new SFBilinearLocal2D(4);
			sfb2d1.asignElement(e);
			sfb2d2.asignElement(e);
			sfb2d3.asignElement(e);
			sfb2d4.asignElement(e);
			mapVS.put(vl.get(0), sfb2d1);
			mapVS.put(vl.get(1), sfb2d2);
			mapVS.put(vl.get(2), sfb2d3);
			mapVS.put(vl.get(3), sfb2d4);			
		} else if(vl.size() == 5) {
			//TODO
			System.out.println("ERROR: getTransformLinear2DShapeFunction() vl.size() == "+vl.size());
		} else {
			System.out.println("ERROR: getTransformLinear2DShapeFunction() vl.size() == "+vl.size());
		}
		return mapVS;
	}	
	
	/**
	 * 利用单元上的形函数进行坐标变换
	 * @param e
	 * @return
	 */
	public Map<Vertex,ShapeFunction> getTransformShapeFunctionByElement(Element e) {
		Map<Vertex, ShapeFunction> mapVS = new LinkedHashMap<Vertex, ShapeFunction>();
		for(int i=1;i<=e.nodes.size();i++) {
			Vertex v = new Vertex(e.nodes.at(i).dim());
			v.set(i,e.nodes.at(i));
			//TODO 用结点上的第一个自由度进行坐标
			mapVS.put(v, e.getDOFList(i).at(1).getShapeFunction());
		}
		return mapVS;
	}
	public List<FunctionDerivable> getTransformFunction(Map<Vertex,ShapeFunction> mapVS) {
		List<FunctionDerivable> r = new LinkedList<FunctionDerivable>();
		int dim = fromVarNames.size();
		for(int i=1;i<=dim;i++) {
			FunctionDerivable a = new FConstant(0.0);
			for(Entry<Vertex,ShapeFunction> e : mapVS.entrySet()) {
				Vertex v = e.getKey();
				ShapeFunction shapeFun = e.getValue();
				FunctionDerivable fm = FOBasicDerivable.Mult(
							new FConstant(v.coord(i)), shapeFun
						);
				a = FOBasicDerivable.Plus(a, fm);
			}
			//???? a.setVarNames(toVarNames);
			r.add(a);
		}
		return r;	
	}	
	
	public void setTransformFunction(List<FunctionDerivable> trans) {
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
	
	public FunctionDerivable getJacobian1D() {
		FunctionDerivable[] funs = new FunctionDerivable[2];
		int index = 0;
		if(transFuns.size() != 2)
			return null;
		for(FunctionDerivable transFun : transFuns) {
			for(String var : toVarNames) {
				DerivativeIndicator di = new DerivativeIndicator();
				di.set(var, 1);
				funs[index++] = transFun.derivative(di);
			}
		}
		/*  
		 *  det(Jac)= Sqrt(0^2 + 1^2)
		 */
		return FOBasicDerivable.Sqrt(FOBasicDerivable.Plus(
				FOBasicDerivable.Mult(funs[0], funs[0]),
				FOBasicDerivable.Mult(funs[1], funs[1]))
				);
	}
	
	public FunctionDerivable getJacobian2D() {
		FunctionDerivable[] funs = new FunctionDerivable[4];
		int index = 0;
		if(transFuns.size() != 2)
			return null;
		for(FunctionDerivable transFun : transFuns) {
			for(String var : toVarNames) {
				DerivativeIndicator di = new DerivativeIndicator();
				di.set(var, 1);
				funs[index++] = transFun.derivative(di);
			}
		}
		/*
		 * 要求结点编号为逆时针  
		 *            |0 1|
		 *  det(Jac)= |2 3| = 0*3-1*2
		 */
		return FOBasicDerivable.Minus(
				FOBasicDerivable.Mult(funs[0], funs[3]),
				FOBasicDerivable.Mult(funs[1], funs[2]));
	}
}
