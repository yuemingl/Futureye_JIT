package edu.uta.futureye.function.basic;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.util.Tools;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.core.TriAreaCoord;
import edu.uta.futureye.core.geometry.Point;
import edu.uta.futureye.core.intf.CoordTrans;
import edu.uta.futureye.function.MultiVarFunc;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.lib.assembler.AssembleParam;
import edu.uta.futureye.util.FutureyeException;
import edu.uta.futureye.util.Utils;

/**
 * <blockquote><pre>
 * Vector to Function
 * Evaluate function values based on vector indices in Variable v
 * 
 * 2011/6/27
 * + Function _d(String varName)
 * 
 * 2011/?/?
 * + extends evaluation ability on the inner area of an element by interpolation 
 * 
 * 2011/10/17
 * + defaultFunction
 * </blockquote></pre>
 * 
 * @author liuyueming
 *
 */
public class Vector2MathFunc extends MultiVarFunc {
	Vector vec = null;
	Mesh mesh = null;
	
	//当自变量超出mesh范围时，使用该函数求值
	MathFunc defaultFunction = null;
	
	boolean enableCache = false;
	Map<Integer, Double> cachedValueMap = null;
	
	/**
	 * 构造一个向量的函数封装，求值时自变量需要指定索引值（=向量的索引），
	 * 自变量如果有坐标值，坐标值将被忽略。
	 * 
	 * @param u
	 */
	public Vector2MathFunc(Vector u) {
		this.vec = u;
	}
	
	/**
	 * <blockquote><pre>
	 * 该构造方法支持导数计算和函数插值，当需要插值计算时，必须指定自变量的坐标值
	 * 
	 * 注意，如果仅需要计算索引值对应的函数值，请使用构造函数public Vector2Function(Vector u)，速度会更快。
	 * 当以下情况发生时，会进行插值计算（假设自变量为Variable v）：
	 * （1）自变量只有坐标值，没有索引值（v.getIndex()==0）
	 * （2）自变量的索引值是其他网格的索引值，不同于当初构造该对象的网格（此时坐标值必须同时指定），
	 *     当使用该构造函数生成的对象作为方程系数进行“整体合成”时，特别是不同网格的插值计算，
	 *     速度会相当慢，建议使用该构造函数构造的对象先插值计算出一个新网格上的Vector，再调用
	 *     构造函数public Vector2Function(Vector u)，构造一个与方程所使用的网格一致的函数对象。
	 * </blockquote></pre>
	 * @param u
	 * @param mesh
	 * @param varName
	 * @param aryVarNames
	 */
	public Vector2MathFunc(Vector u, Mesh mesh,
			String varName, String ...aryVarNames) {
		this.vec = u;
		this.mesh = mesh;
		if(vec.getDim() != mesh.getNodeList().size())
			throw new FutureyeException("u.getDim() != mesh.getNodeList().size()");
		varNames = new String[1+aryVarNames.length];
		varNames[0] = varName;
		for(int i=0; i<aryVarNames.length; i++)
			varNames[i+1] = aryVarNames[i];
	}
	
	public Vector2MathFunc(Vector u, Mesh mesh,
			String []aryVarNames) {
		this.vec = u;
		this.mesh = mesh;
		if(vec.getDim() != mesh.getNodeList().size())
			throw new FutureyeException("u.getDim() != mesh.getNodeList().size()");
		varNames = aryVarNames;
	}
	
	public Vector2MathFunc(Vector u, Mesh mesh,
			List<String> varNames) {
		this.vec = u;
		this.mesh = mesh;
		this.varNames = varNames.toArray(new String[0]);
	}	
	
	/**
	 * 当自变量超出mesh范围时，使用该函数求值
	 * @param fun
	 * @return
	 */
	public Vector2MathFunc setDefaultFunction(MathFunc fun) {
		this.defaultFunction = fun;
		return this;
	}
	public MathFunc getDefaultFunction() {
		return this.defaultFunction;
	}
	
	/**
	 * 缓存插值出来的函数值，以加快相同坐标值重复计算的速度
	 * 注意：使用前最好先清除缓存
	 * @param b
	 */
	public void cacheInterpolateValue(boolean b) {
		this.enableCache = b;
		if(b && cachedValueMap == null)
			cachedValueMap = new HashMap<Integer, Double>();
		else if(!b)
			cachedValueMap = null;
	}
	/**
	 * 清除插值出来的函数值缓存
	 */
	public void clearInterpolateValueCache() {
		if(cachedValueMap != null)
			cachedValueMap.clear();
	}
	public Mesh getMesh() {
		return this.mesh;
	}
	public Vector getVector() {
		return this.vec;
	}
	
	@Override
	public double apply(Variable v) {
		int index = v.getIndex();
		int nDim = vec.getDim();
		if(mesh == null) { //完全依靠index来求值
			if(index > 0 && index <= nDim)
				return vec.get(index);//注：下标错位会造成结果出现随机混乱
			else if(index == 0) {
				throw new FutureyeException("Error: index(=0), please specify index of Variable!");
			} else
				throw new FutureyeException(String.format("Error: index(=%d) should between [1,%d]\n"+
						"try to use constructor Vector2Function(Vector u, Mesh mesh, String varName, String ...aryVarNames)\n"+
						"to extends evaluation ability on the inner area of an element by interpolation.",
						index,nDim));
		} else { //如果同时指定mesh，需要判断是否需要插值
			boolean needInterpolation = false;
			double[] coord = new double[varNames.length];
			for(int i=0;i<varNames.length;i++) {
				coord[i] = v.get(varNames[i]);
			}
			
			if(index < 0)
				throw new FutureyeException(String.format("Error: index(=%d) <0!",index));
			else if(index==0 || index > nDim)
				needInterpolation = true;
			else {
				//判断是否需要插值，以便节省计算时间
				Node node = mesh.getNodeList().at(index);
				Node node2 = new Node().set(0,coord);
				//判断v的坐标与mesh.getNodeList().at(index)是否一致，不一致需要插值
				if( ! node.coordEquals(node2) )
					needInterpolation = true;
			}
			
			if(needInterpolation) {
				if(enableCache && index > 0) {
					Double cacheValue = cachedValueMap.get(index);
					if(cacheValue != null) return cacheValue;
				}

				//先找到包含v的坐标的单元，在单元内进行插值
				Element e = mesh.getElementByCoord(coord);
				if(e == null) {
					if(this.defaultFunction == null)
						throw new FutureyeException(
								"Can't interplate "+v.toString()+" on the current mesh! \n"+
								"Try to use setDefaultFunction(Function fun) to extend the evaluation outside of the mesh domain.");
					else
						return this.defaultFunction.apply(v);
				}
				double[] f = new double[e.nodes.size()];
				for(int i=1;i<=e.nodes.size();i++) {
					f[i-1] = vec.get(e.nodes.at(i).globalIndex);
				}
				//二维四边形单元
				if(e.vertices().size() == 4 && coord.length==2) {
					double[] coef = Utils.computeBilinearFunctionCoef(e.nodes.toArray(new Point[0]), f);
					//f(x,y) = a1 + a2*x + a3*y + a4*x*y
					double x = coord[0];
					double y = coord[1];
					double interpValue = coef[0] + coef[1]*x + coef[2]*y + coef[3]*x*y;
					if(enableCache && index > 0) {
						cachedValueMap.put(index, interpValue);
					}
					return interpValue;
				}
				throw new FutureyeException("Error: Unsported element type:"+e.toString());
			} else {
				return vec.get(index);
			}
		}
	}
	

	@Override
	public double apply(AssembleParam ap, double... args) {
		
		Element e = ap.element;
		int i1 = e.nodes.at(1).globalIndex;
		int i2 = e.nodes.at(2).globalIndex;
		int i3 = e.nodes.at(3).globalIndex;
		double v1 = this.vec.get(i1);
		double v2 = this.vec.get(i2);
		double v3 = this.vec.get(i3);
		int startIdx = 12;
		return v1*args[startIdx] + v2*args[startIdx+1] + v3*args[startIdx+2];
		
	}

	@Override
	public double apply(double... args) {
		return apply(null, args);
	}

	
	@Override
	public MathFunc diff(String varName) {
		if(mesh == null) 
			throw new FutureyeException(
					"Please use constructor Vector2Function(Vector u, Mesh mesh, String varName, String ...aryVarNames)");
		Vector vd = Tools.computeDerivative(mesh, vec, varName);
		MathFunc fd = new Vector2MathFunc(vd,mesh,this.varNames);
		return fd;
	}
	
	@Override
	public String toString() {
		return "Vector2Function:"+this.vec.toString();
	}

}
