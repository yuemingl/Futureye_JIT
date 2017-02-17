package edu.uta.futureye.util;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import edu.uta.futureye.algebra.SpaceVector;
import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.core.CoordinateTransform;
import edu.uta.futureye.core.DOF;
import edu.uta.futureye.core.DOFOrder;
import edu.uta.futureye.core.Edge;
import edu.uta.futureye.core.EdgeLocal;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Face;
import edu.uta.futureye.core.FaceLocal;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.core.NodeRefined;
import edu.uta.futureye.core.NodeType;
import edu.uta.futureye.core.Vertex;
import edu.uta.futureye.core.Volume;
import edu.uta.futureye.core.geometry.GeoEntity;
import edu.uta.futureye.core.geometry.GeoEntity2D;
import edu.uta.futureye.core.geometry.GeoEntity3D;
import edu.uta.futureye.core.geometry.Point;
import edu.uta.futureye.core.intf.FiniteElement;
import edu.uta.futureye.function.MultiVarFunc;
import edu.uta.futureye.function.FMath;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.basic.SpaceVectorFunction;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.function.intf.ScalarShapeFunction;
import edu.uta.futureye.function.intf.VectorMathFunc;
import edu.uta.futureye.util.container.DOFList;
import edu.uta.futureye.util.container.ElementList;
import edu.uta.futureye.util.container.NodeList;
import edu.uta.futureye.util.container.ObjIndex;
import edu.uta.futureye.util.container.ObjList;
import edu.uta.futureye.util.container.VertexList;

public class Utils {
	
	public static List<String> mergeList(List<String> a, List<String> b) {
		Set<String> set = new LinkedHashSet<String>();
		if(a != null)
			set.addAll(a);
		if(b != null)
			set.addAll(b);
		List<String> rlt = new LinkedList<String>();
		rlt.addAll(set);
		Collections.sort(rlt, new Comparator<String>() {
			@Override
			public int compare(String o1, String o2) {
				return o1.compareTo(o2);
			}
			
		});
		return rlt;
	}
	
	public static Map<String, Integer> getIndexMap(List<String> list) {
		Map<String, Integer> map = new HashMap<String, Integer>();
		for(int i=0; i<list.size(); i++) {
			map.put(list.get(i), i);
		}
		return map;
	}
	
	public static boolean isMapContain(Map<String, Integer> container, Map<String, Integer> target) {
		if(target == null) return false;
		if(container.size() < target.size()) return false;
		for(Entry<String, Integer> e : target.entrySet()) {
			if(container.get(e.getKey()) != e.getValue())
				return false;
		}
		return true;
	}
	
	public static boolean isListEqualIgnoreOrder(List<String> l1, List<String> l2) {
		if(l1.size() != l2.size()) return false;
		List<String> t1 = new ArrayList<String>();
		t1.addAll(l1);
		List<String> t2 = new ArrayList<String>();
		t2.addAll(l2);
		Collections.sort(t1, new Comparator<String>() {
			@Override
			public int compare(String o1, String o2) {
				return o1.compareTo(o2);
			}
		});
		Collections.sort(t2, new Comparator<String>() {
			@Override
			public int compare(String o1, String o2) {
				return o1.compareTo(o2);
			}
		});
		
		return t1.equals(t2); //
	}
	
	/**
	 * TODO 可以修改为2D,3D都可以计算
	 */
	public static double computeAngle2D(Point a1,Point a2,Point b1,Point b2) {
		Vector v1 = new SpaceVector(2);
		v1.set(1, a2.coord(1)-a1.coord(1));
		v1.set(2, a2.coord(2)-a1.coord(2));
		
		Vector v2 = new SpaceVector(2);
		v2.set(1, b2.coord(1)-b1.coord(1));
		v2.set(2, b2.coord(2)-b1.coord(2));
		
		//if(v1.norm2()<Constant.eps || v2.norm2()<Constant.eps), Math.acos(v) will be NaN

		double v = v1.dot(v2)/(v1.norm2()*v2.norm2());
		if(v > 1.0) v = 1.0;
		else if(v < -1.0) v = -1.0;
		return Math.acos(v);
		
	}

	public static double computeAngle(SpaceVector v1,SpaceVector v2) {
		double v = v1.dot(v2)/(v1.norm2()*v2.norm2());
		if(v > 1.0) v = 1.0;
		else if(v < -1.0) v = -1.0;
		return Math.acos(v);
	}
	
	/**
	 * y1=f(x1)
	 * y2=f(x2)
	 * ==>
	 * return y=f(x)
	 * @param x1
	 * @param x2
	 * @param x
	 * @param y1
	 * @param y2
	 * @return
	 */
	public static double linearInterpolate(double x1, double x2, double x, 
			double y1, double y2) {
		double k = (y2 -y1)/(x2 -x1);
		double r = k*(x-x1)+y1;
		return r;
	}
	
	
	/**
	 * y1=f(p1)
	 * y2=f(p2)
	 * ==>
	 * return y=f(p)
	 * TODO 返回p?
	 * @param p1
	 * @param p2
	 * @param p
	 * @param y1
	 * @param y2
	 * @return
	 */
	public static double linearInterpolate(Point p1, Point p2, Point p, 
			double y1, double y2) {
		double x1,x2,x;
		for(int i=1;i<=p1.dim();i++) {
			x1 = p1.coord(i);
			x2 = p2.coord(i);
			x = p.coord(i);
			if(Math.abs(x1-x2)>Constant.meshEps) {
				return linearInterpolate(x1,x2,x,y1,y2);
			}
		}
		return 0.0;
	}
	
	public static Double getRefCoord(double a,double b,double x) {
		double t = 0.0;
		if(Math.abs(a-b)<Constant.eps) {
			return null;
		} else if(a > b) {
			t=b; b=a; a=t;
		}
		double rlt = 2.0*(x-a)/(b-a)-1.0;
		return rlt;
	}
	
	public static double quadraticInterpolate(Point p1, Point p2, Point p3, Point p, 
			double y1, double y2, double y3) {
//		SFQuadraticLocal1D sp1 = new SFQuadraticLocal1D(1);
//		SFQuadraticLocal1D sp2 = new SFQuadraticLocal1D(2);
//		SFQuadraticLocal1D sp3 = new SFQuadraticLocal1D(3);
		Double x = 0.0;
		int i;
		for(i=1;i<=p1.dim();i++) {
			x = getRefCoord(p1.coord(i),p3.coord(i),p.coord(i));
			if(x != null) break;
		}
		double rlt;
		if(x < 0) {
			if(p1.coord(i) < p3.coord(i))
				rlt = linearInterpolate(p1.coord(i),p2.coord(i),p.coord(i),y1,y2);
			else
				rlt = linearInterpolate(p3.coord(i),p2.coord(i),p.coord(i),y3,y2);
		} else {
			if(p1.coord(i) < p3.coord(i))
				rlt = linearInterpolate(p2.coord(i),p3.coord(i),p.coord(i),y2,y3);
			else
				rlt = linearInterpolate(p2.coord(i),p1.coord(i),p.coord(i),y2,y1);
		}
//		double x1 = p1.coord(i);
//		double x2 = p3.coord(i);
//		double ytmp = 0.0;
//		if(x1 > x2) {
//			ytmp=y1; y1=y3; y3=ytmp; 
//		}
//		Variable var = new Variable();
//		var.set("r", x);
//		double rlt = y1*sp1.value(var)+y2*sp2.value(var)+y3*sp3.value(var);
		return rlt;
	}	
	
	/**
	 * Gauss smooth
	 * 
	 * @param mesh
	 * @param u
	 * @param neighborBand 1: one neighbor; 2 layers of neighbor; 3...
	 * @param weight weight*curNode + (1-weight)*neighborNode
	 * @return
	 */
	public static <Vec extends Vector> Vec gaussSmooth(Mesh mesh, Vec u, int neighborBand, double weight) {
		NodeList list = mesh.getNodeList();
		int nNode = list.size();
		if(nNode != u.getDim()) {
			throw new FutureyeException("Node number of mesh != length of vector u: "+
					"nNode="+nNode+"  dim(u)="+u.getDim());
		}
		
	    @SuppressWarnings("unchecked")
		Vec su = (Vec)u.copy();
	    su.setAll(0.0);
	    
	    if(neighborBand == 1) {
		    for(int i=1;i<=nNode;i++) {
		    	NodeList nbList = new NodeList();
		    	Node node = list.at(i);
		    	//TODO 自适应网格节点需要注意
		    	if(node instanceof NodeRefined) {
		    		if(((NodeRefined) node).isHangingNode()) {
		    			NodeList cns = ((NodeRefined) node).constrainNodes;
		    			for(int k=1;k<=cns.size();k++) {
		    				nbList.addAll(cns.at(k).neighbors);
		    			}
		    		}
		    	}
		    	nbList.addAll(node.neighbors);
		    	if(nbList==null || nbList.size() == 0) {
					throw new FutureyeException("No beighbors of Node "+node.globalIndex+
							", call mesh.computeNeiborNode() first!");
		    	}
		    	double nbV = 0.0;
		    	for(int j=1;j<=nbList.size();j++) {
		    		Node nbNode = nbList.at(j);
		    		nbV += u.get(nbNode.globalIndex);
		    	}
		    	double rlt = weight*u.get(node.globalIndex)+(1-weight)*(nbV/nbList.size());
		    	su.set(node.globalIndex, rlt);
		    }
	    } else {
	    	Map<Node, Integer> map = new HashMap<Node, Integer>();
		    for(int i=1;i<=nNode;i++) {
		    	Node node = list.at(i);
		    	map.clear();
		    	map.put(node, 0);

		    	NodeList nbList = node.neighbors;
		    	if(nbList==null || nbList.size() == 0) {
					throw new FutureyeException("No beighbors of Node "+node.globalIndex+
							", call mesh.computeNeiborNode() first!");
		    	}
		    	for(int j=1;j<=nbList.size();j++) {
		    		Node nbNode = nbList.at(j);
		    		map.put(nbNode, 1);
		    		
		    		NodeList nb2List = nbNode.neighbors;
		    		for(int k=1;k<=nb2List.size();k++) {
		    			Node nb2Node = nb2List.at(k);
		    			if(map.get(nb2Node) == null) 
		    				map.put(nb2Node, 2);
		    		}
		    	}
		    	int nbNodeNum = 0;
		    	int nb2NodeNum = 0;
		    	double nbValue = 0.0;
		    	double nb2Value = 0.0;
		    	for(Entry<Node,Integer> e : map.entrySet()) {
		    		Node n = e.getKey();
		    		if(e.getValue() == 1) {
		    			nbNodeNum++;
		    			nbValue += u.get(n.globalIndex);
		    		} else if(e.getValue() == 2) {
		    			nb2NodeNum++;
		    			nb2Value += u.get(n.globalIndex);
		    		}
		    	}
		    	double rlt = weight*u.get(node.globalIndex)+
		    	(1-weight)*weight*(nbValue/nbNodeNum)+
		    	(1-weight)*(1-weight)*(nb2Value/nb2NodeNum);
		    	su.set(node.globalIndex, rlt);
		    }
	    }
	    return su;
	}
	
	/**
	 * 每个结点的值取为周围结点的最大值
	 * @param mesh
	 * @param u
	 * @return
	 */
	public static Vector gaussMax(Mesh mesh, Vector u) {
		NodeList list = mesh.getNodeList();
		int nNode = list.size();
		if(nNode != u.getDim()) {
			throw new FutureyeException("Node number of mesh != length of vector u: "+
					"nNode="+nNode+"  dim(u)="+u.getDim());
		}
		
	    Vector su = u.copy();
	    su.setAll(0.0);
	    
	    for(int i=1;i<=nNode;i++) {
	    	NodeList nbList = new NodeList();
	    	Node node = list.at(i);
	    	//TODO 自适应网格节点需要注意
	    	if(node instanceof NodeRefined) {
	    		if(((NodeRefined) node).isHangingNode()) {
	    			NodeList cns = ((NodeRefined) node).constrainNodes;
	    			for(int k=1;k<=cns.size();k++) {
	    				nbList.addAll(cns.at(k).neighbors);
	    			}
	    		}
	    	}
	    	nbList.addAll(node.neighbors);
	    	if(nbList.size() == 0) {
				throw new FutureyeException("No beighbors of Node "+node.globalIndex+
						", call mesh.computeNeiborNode() first!");
	    	}
	    	double nbV = Double.MIN_VALUE;
	    	for(int j=1;j<=nbList.size();j++) {
	    		Node nbNode = nbList.at(j);
	    		if(nbV < u.get(nbNode.globalIndex))
	    			nbV = u.get(nbNode.globalIndex);
	    	}
	    	double rlt = nbV;
	    	su.set(node.globalIndex, rlt);
	    }
	    return su;
	}
	
	public static VectorMathFunc interpolateOnElement(VectorMathFunc fun, Element e) {
		SpaceVectorFunction rlt = new SpaceVectorFunction(fun.getDim());
		for(int i=1;i<=fun.getDim();i++) {
			rlt.set(i, interpolateOnElement(fun.get(i),e));
		}
		return rlt;
	}
	public static MathFunc interpolateOnElement(MathFunc fun, Element e) {
		if(fun == null) throw new FutureyeException("fun should not be null");
		if(fun.isConstant())
			return fun;
		MathFunc rlt = new FC(0.0);
//		int nNode = e.nodes.size();
//		for(int i=1;i<=nNode;i++) {
//			DOFList dofListI = e.getNodeDOFList(i);
//			for(int k=1;k<=dofListI.size();k++) {
//				DOF dofI = dofListI.at(k);
//				Variable var = Variable.createFrom(fun, (Node)dofI.getOwner(), dofI.getGlobalNumber());
//				Function PValue = new FConstant(fun.value(var));
//				rlt = FOBasic.Plus(rlt, FOBasic.Mult(PValue, dofI.getSSF()));
//			}
//		}
		if(e.dim() == 1) {
			CoordinateTransform trans = new CoordinateTransform(1);
			Map<Vertex,ScalarShapeFunction> transSF = trans.getTransformLinear1DShapeFunction(e);
			for(Entry<Vertex,ScalarShapeFunction> entry : transSF.entrySet()) {
				Point p = entry.getKey();
				ScalarShapeFunction sf = entry.getValue();
				int index = 0;
				//当fun是向量值函数时，函数的值靠下标获得，应该获取结点上自由度的全局编号
				int localNodeIndex = e.getLocalIndex(e.getNode(p));
				DOFList DOFs = e.getNodeDOFList(localNodeIndex);
				if(DOFs != null)
					index = DOFs.at(1).getGlobalIndex();
				Variable var = Variable.createFrom(fun, p, index);
				var.setElement(e);
				MathFunc PValue = new FC(fun.apply(var));
				rlt = rlt.A(PValue.M(sf));			
			}
			
		} else if(e.dim() == 2) {
			CoordinateTransform trans = new CoordinateTransform(2);
			Map<Vertex,ScalarShapeFunction> transSF = trans.getTransformLinear2DShapeFunction(e);
			for(Entry<Vertex,ScalarShapeFunction> entry : transSF.entrySet()) {
				Point p = entry.getKey();
				ScalarShapeFunction sf = entry.getValue();
				int index = 0;
				//当fun是向量值函数时，函数的值靠下标获得，应该获取结点上自由度的全局编号
				int localNodeIndex = e.getLocalIndex(e.getNode(p));
				DOFList DOFs = e.getNodeDOFList(localNodeIndex);
				if(DOFs != null)
					index = DOFs.at(1).getGlobalIndex();
				Variable var = Variable.createFrom(fun, p, index);
				var.setElement(e);
				MathFunc PValue = new FC(fun.apply(var));
				rlt = rlt.A(PValue.M(sf));	
			}
		} else if(e.dim()==3) {
			//TODO 有没有更简洁的办法？
			CoordinateTransform trans = new CoordinateTransform(3);
			Map<Vertex,ScalarShapeFunction> transSF = trans.getTransformLinear3DShapeFunction(e);
			for(Entry<Vertex,ScalarShapeFunction> entry : transSF.entrySet()) {
				Point p = entry.getKey();
				ScalarShapeFunction sf = entry.getValue();
				int index = 0;
				//当fun是向量值函数时，函数的值靠下标获得，应该获取结点上自由度的全局编号
				int localNodeIndex = e.getLocalIndex(e.getNode(p));
				DOFList DOFs = e.getNodeDOFList(localNodeIndex);
				if(DOFs != null)
					index = DOFs.at(1).getGlobalIndex();
				Variable var = Variable.createFrom(fun, p, index);
				var.setElement(e);
				MathFunc PValue = new FC(fun.apply(var));
				rlt = rlt.A(PValue.M(sf));			
			}
		}
		return rlt;
	}
	
	//只是用于三角形线性元
	public static Map<String, MathFunc> getFunctionComposeMap(Element e) {
		final Element fe =e;
		Map<String, MathFunc> fInners = new HashMap<String, MathFunc>();
		final List<String> varNamesInner = new LinkedList<String>();
		varNamesInner.add("r");
		varNamesInner.add("s");
		varNamesInner.add("t");
		fInners.put("x", new MultiVarFunc("x",varNamesInner) {	
//			@Override
//			public double apply(Variable v) {
//				double rlt = 0.0;
//				for(int i=1;i<=fe.nodes.size();i++)
//					rlt += fe.nodes.at(i).coord(1)*v.get(varNamesInner.get(i-1));
//				return rlt;
//			}
			
			@Override
			public double apply(double... args) {
				double rlt = 0.0;
				for(int i=1;i<=fe.nodes.size();i++)
					rlt += fe.nodes.at(i).coord(1)*args[this.argIdx[i-1]];
				return rlt;
			}
		});
		fInners.put("y", new MultiVarFunc("y", varNamesInner) {	
//			@Override
//			public double apply(Variable v) {
//				double rlt = 0.0;
//				for(int i=1;i<=fe.nodes.size();i++)
//					rlt += fe.nodes.at(i).coord(2)*v.get(varNamesInner.get(i-1));
//				return rlt;
//			}
			
			@Override
			public double apply(double... args) {
				double rlt = 0.0;
				for(int i=1;i<=fe.nodes.size();i++)
					rlt += fe.nodes.at(i).coord(1)*args[this.argIdx[i-1]];
				return rlt;
			}
		});		
		return fInners;
	}
	
	/**
	 * 判断点p是否在p1,p2直线段中间，包括端点
	 * @param p1
	 * @param p2
	 * @param p
	 * @return
	 */
	public static boolean isPointOnLineSegment(Point p1, Point p2, Point p) {
		if(p1.coordEquals(p) || p2.coordEquals(p))
			return true;
		double rlt = computeAngle2D(p, p1, p, p2);
		if(Math.abs(rlt-Math.PI) < Constant.meshEps) {
			return true;
		}
		return false;
	}
	
	/**
	 * 判断点p是否在p1,p2直线段中间，不包括端点
	 * @param p1
	 * @param p2
	 * @param p
	 * @return
	 */
	public static boolean isPointOnLineSegmentNoEndingPoint(Point p1, Point p2, Point p) {
		if(p1.coordEquals(p) || p2.coordEquals(p))
			return false;
		double rlt = computeAngle2D(p, p1, p, p2);
		if(Math.abs(rlt-Math.PI) < Constant.meshEps) {
			return true;
		}
		return false;
	}
	
	/**
	 * 判断点p是否在p1,p2直线上，不一定在两点之间
	 * @param p1
	 * @param p2
	 * @param p
	 * @return
	 */
	public static boolean isPointOnLine(Point p1, Point p2, Point p) {
		if(p1.coordEquals(p) || p2.coordEquals(p))
			return true;
		double rlt = computeAngle2D(p, p1, p, p2);
		if(Math.abs(rlt-Math.PI) < Constant.meshEps ||
				Math.abs(rlt) < Constant.meshEps) {
			return true;
		}
		return false;
	}
	
	/**
	 * 判断线段[a1,a2] 与 [b1,b2]是否有公共线段
	 * @param a1
	 * @param a2
	 * @param b1
	 * @param b2
	 * @return
	 */
	public static boolean isLineOverlap(Point a1,Point a2,Point b1,Point b2) {
		if(a1.coordEquals(b1) && a2.coordEquals(b2))
			return true;
		if(a2.coordEquals(b1) && a1.coordEquals(b2))
			return true;
		if(isPointOnLineSegmentNoEndingPoint(a1,a2,b1) ||isPointOnLineSegmentNoEndingPoint(a1,a2,b2) ||
				isPointOnLineSegmentNoEndingPoint(b1,b2,a1) || isPointOnLineSegmentNoEndingPoint(b1,b2,a2) )
			return true;
		return false;
	}
	
	/**
	 * 计算两个点之间的距离
	 * 
	 * TODO 考虑是否在Point上建立代数运算方法？
	 * 
	 * @param p1
	 * @param p2
	 * @return
	 */
	public static double computeLength(Point p1,Point p2) {
		int dim = p1.dim();
		double len = 0.0;
		for(int i=1;i<=dim;i++) {
			double t = (p1.coord(i)-p2.coord(i));
			len += t*t;
		}
		return Math.sqrt(len);
	}
	
	/**
	 * Return a unit outer norm vector of a 2D edge (a,b), where
	 * n \dot (a-b) = 0
	 * |b| = 1
	 * @param edge
	 * @return n
	 */
	public static Vector getNormVector(Point a,Point b) {
		SpaceVector v1 = new SpaceVector(a.coords());
		SpaceVector v2 = new SpaceVector(b.coords());
		Vector v2mv1 = FMath.axpy(-1.0, v1, v2);//v2 - v1
		double len = v2mv1.norm2();
		Vector rlt =  FMath.ax(1.0/len, v2mv1);
		return new SpaceVector(rlt.get(2),-rlt.get(1));//外法向
	}
	
	/**
	 * 计算三角形面积（二维点或者三维点）
	 * @param vertices
	 * @return
	 */
	public static double getTriangleArea(ObjList<Vertex> vertices) {
		double area = 0.0;
		if(vertices.size() == 3 && vertices.at(1).dim() == 2) {
			double x1 = vertices.at(1).coord(1) , y1 =  vertices.at(1).coord(2) ;
			double x2 = vertices.at(2).coord(1) , y2 =  vertices.at(2).coord(2) ;
			double x3 = vertices.at(3).coord(1) , y3 =  vertices.at(3).coord(2) ;
			area = ( (x2*y3 - x3*y2) - (x1*y3 - x3*y1) + (x1*y2 - x2*y1) ) / 2.0;
		} else if(vertices.size() == 3 && vertices.at(1).dim() == 3) {
			double x1 = vertices.at(1).coord(1) , y1 =  vertices.at(1).coord(2) , z1 =  vertices.at(1).coord(3) ;
			double x2 = vertices.at(2).coord(1) , y2 =  vertices.at(2).coord(2) , z2 =  vertices.at(2).coord(3) ;
			double x3 = vertices.at(3).coord(1) , y3 =  vertices.at(3).coord(2) , z3 =  vertices.at(3).coord(3) ;
			SpaceVector v1 = new SpaceVector(3);
			v1.set(1, x2-x1);
			v1.set(2, y2-y1);
			v1.set(3, z2-z1);
			SpaceVector v2 = new SpaceVector(3);
			v2.set(1, x3-x1);
			v2.set(2, y3-y1);
			v2.set(3, z3-z1);
			area = v1.crossProduct(v2).norm2() / 2.0;
			if(Math.abs(area)<Constant.eps) throw new FutureyeException();
		} else {
			throw new FutureyeException();
		}
		if(Math.abs(area)<Constant.eps) throw new FutureyeException();
		return area;
	}
	
	/**
	 * 计算二维四边形面积
	 * @param vertices
	 * @return
	 */
	public static double getRectangleArea(ObjList<Vertex> vertices) {
		double area = 0.0;
		if(vertices.size() == 4) {
			area += getTriangleArea(vertices.subList(1, 3));
			area += getTriangleArea(vertices.subList(new ObjIndex(1,3,4)));
		} else {
			throw new FutureyeException("Number of vertices != 4");
		}
		return area;
	}
	
	/**
	 * 多边形（二维）面积，转化为计算多个三角形面积
	 * @param vertices
	 * @return
	 */
	public static double getPolygonArea(ObjList<Vertex> vertices) {
		double area = 0.0;
		for(int i=3;i<=vertices.size();i++) {
			area += getTriangleArea(
					vertices.subList(new ObjIndex(1,i-1,i))
					);
		}
		return area;
	}
	
	/**
	 * 计算球面三角形面积
	 * @param r
	 * @param o center
	 * @param a
	 * @param b
	 * @param c
	 * @return
	 */
	public static double getSphereTriangleArea(double r,Point o, Point a, Point b, Point c) {
		
		SpaceVector oa = new SpaceVector(
				a.coord(1)-o.coord(1),
				a.coord(2)-o.coord(2),
				a.coord(3)-o.coord(3)
				);
		SpaceVector ob = new SpaceVector(
				b.coord(1)-o.coord(1),
				b.coord(2)-o.coord(2),
				b.coord(3)-o.coord(3)
				);		
		SpaceVector oc = new SpaceVector(
				c.coord(1)-o.coord(1),
				c.coord(2)-o.coord(2),
				c.coord(3)-o.coord(3)
				);
		
		SpaceVector v1,v2,v3;
		v1 = oa.crossProduct(ob);
		v2 = ob.crossProduct(oc);
		v3 = oc.crossProduct(oa);
		double angle1 = Math.PI-computeAngle(v1,v2);
		double angle2 = Math.PI-computeAngle(v2,v3);
		double angle3 = Math.PI-computeAngle(v3,v1);
//Ouput Matlab code to visualize parameters
//		System.out.print("plot3(");
//		System.out.print("["+o.coord(1)+" "+a.coord(1)+" "+b.coord(1)+" "+c.coord(1)+" "+a.coord(1)+"],");
//		System.out.print("["+o.coord(2)+" "+a.coord(2)+" "+b.coord(2)+" "+c.coord(2)+" "+a.coord(2)+"],");
//		System.out.print("["+o.coord(3)+" "+a.coord(3)+" "+b.coord(3)+" "+c.coord(3)+" "+a.coord(3)+"])");
//		System.out.print("\nhold on\n");
//		System.out.print("plot3(");
//		System.out.print("["+b.coord(1)+" "+o.coord(1)+" "+c.coord(1)+"],");
//		System.out.print("["+b.coord(2)+" "+o.coord(2)+" "+c.coord(2)+"],");
//		System.out.print("["+b.coord(3)+" "+o.coord(3)+" "+c.coord(3)+"])");
//		System.out.println("\nhold on\n");
		
		double area = (angle1 + angle2 + angle3 - Math.PI)*r*r;

		return area;
	}
	

	public static double getTetrahedronVolume(ObjList<Vertex> vertices) {
		double x1,x2,x3,x4;
		double y1,y2,y3,y4;
		double z1,z2,z3,z4;
		
		x1 = vertices.at(1).coord(1);
		x2 = vertices.at(2).coord(1);
		x3 = vertices.at(3).coord(1);
		x4 = vertices.at(4).coord(1);
		y1 = vertices.at(1).coord(2);
		y2 = vertices.at(2).coord(2);
		y3 = vertices.at(3).coord(2);
		y4 = vertices.at(4).coord(2);
		z1 = vertices.at(1).coord(3);
		z2 = vertices.at(2).coord(3);
		z3 = vertices.at(3).coord(3);
		z4 = vertices.at(4).coord(3);
		
		/*
		      |x2-x1 x3-x1 x4-x1| |1 2 3|
		6*v = |y2-y1 y3-y1 y4-y1|=|4 5 6|=1*(5*9-8*6) + 4*(8*3-2*9) + 7*(2*9-8*3)
		      |z2-z1 z3-z1 z4-z1| |7 8 9|
		*/
		double volume = (x2-x1)*((y3-y1)*(z4-z1)-(y4-y1)*(z3-z1))
			   + (y2-y1)*((x4-x1)*(z3-z1)-(x3-x1)*(z4-z1))
			   + (z2-z1)*((x3-x1)*(y4-y1)-(x4-x1)*(y3-y1));
		volume = volume/6.0;
		return volume;
	}
	
	public static double getHexahedronVolume(ObjList<Vertex> vertices) {
		return 0.0;
		
	}
	
	public static double getPolyhedronVolume(ObjList<Vertex> vertices) {
		return 0.0;
		
	}
	
	/**
	 * f(x,y) = a1 + a2*x + a3*y + a4*x*y
	 * 
	 * @param p: [(x1,y1) (x2,y2) (x3,y3) (x4,y4)]
	 * @param f: [f1 f2 f3 f4]
	 * @return [a1 a2 a3 a4]
	 */
	public static double[] computeBilinearFunctionCoef(
				Point[] p, double[] fValues) {
		int len = p.length;
		if(len != 4) {
			throw new FutureyeException("p.size()="+len+", should be 4.");
		}
		
		
		double[] x = new double[5];//x [1...4]
		double[] y = new double[5];//y [1...4]
		double[] xx= new double[5];//xx[1...4]
		double[] yy= new double[5];//yy[1...4]
		double[] f = new double[5];//f [1...4]
		double[] a = new double[4];//a [0...3]
		for(int i=1;i<=4;i++) {
			x[i] = p[i-1].coord(1);
			y[i] = p[i-1].coord(2);
			xx[i] = p[i-1].coord(1);
			yy[i] = p[i-1].coord(2);
			f[i] = fValues[i-1];
		}

		int r = -1;
		for(int i=2;i<=len;i++) {
			if(x[1] != x[i]) {
				r = i;
				break;
			}
		}
		if(r < 0) {
			throw new FutureyeException("r<0");
		}
		if(r != 2) {
			x[2] = xx[r];
			x[r] = xx[2];
			y[2] = yy[r];
			y[r] = yy[2];
			f[2] = fValues[r-1];
			f[r] = fValues[1];
		}
		
		
		double x21 = x[2]-x[1];
		double x31 = x[3]-x[1];
		double x41 = x[4]-x[1];
		double y21 = y[2]-y[1];
		double y31 = y[3]-y[1];
		double y41 = y[4]-y[1];
		double x1y1 = x[1]*y[1];
		double x2y2 = x[2]*y[2];
		double x3y3 = x[3]*y[3];
		double x4y4 = x[4]*y[4];
		
		/**
		 *  (b11 b12)(a3) = (g1)
		 *  (b21 b22)(a4)   (g2)
		 */
		double b11 = x21*y31 - y21*x31;
		double b12 = x21*(x3y3-x1y1) - (x2y2-x1y1)*x31;
		double b21 = x21*y41 - y21*x41;
		double b22 = x21*(x4y4-x1y1) - (x2y2-x1y1)*x41;
		double g1  = x21*(f[3]-f[1]) - (f[2]-f[1])*x31;
		double g2  = x21*(f[4]-f[1]) - (f[2]-f[1])*x41;
		
		double bbbb = b11*b22-b21*b12;
		a[2] = (b22*g1-b12*g2)/bbbb;
		a[3] = (b11*g2-b21*g1)/bbbb;
		
		a[1] = ((f[2]-f[1]) - (y[2]-y[1])*a[2] - (x2y2-x1y1)*a[3])/(x21);
		a[0] = f[1] - x[1]*a[1] - y[1]*a[2] - x1y1*a[3];
		
//		if(r != 2) {
//			double tmp = a[r-1];
//			a[r-1] = a[1];
//			a[1] = tmp;
//		}
		return a;
	}
	
	public static double determinant(double[][] A) {
		if(A.length == 2) {
			return A[0][0]*A[1][1] - A[0][1]*A[1][0];
		} else if(A.length == 3) {
			return    A[0][0]*(A[1][1]*A[2][2] - A[1][2]*A[2][1]) 
					- A[0][1]*(A[1][0]*A[2][2] - A[1][2]*A[2][0]) 
					+ A[0][2]*(A[0][0]*A[1][1] - A[0][1]*A[0][1]);
		}
		throw new FutureyeException("NOT Supported!");
	}

	/**
	 * 
	 * @param A A[2][2][len] or A[3][3][len]
	 * @return
	 */
	public static double[] determinant(double[][][] A) {
		int len = A[0][0].length;
		double[] rlt = new double[len];
		if(A.length == 2) {
			for(int i=0;i<len;i++) {
					rlt[i] = A[0][0][i]*A[1][1][i] - A[0][1][i]*A[1][0][i];
			}
		} else if(A.length == 3) {
			for(int i=0;i<len;i++) {
				rlt[i] = A[0][0][i]*(A[1][1][i]*A[2][2][i] - A[1][2][i]*A[2][1][i]) 
						- A[0][1][i]*(A[1][0][i]*A[2][2][i] - A[1][2][i]*A[2][0][i]) 
						+ A[0][2][i]*(A[0][0][i]*A[1][1][i] - A[0][1][i]*A[0][1][i]);
			}
		} else 
			throw new FutureyeException("NOT Supported!");
		return rlt;
	}
	
	/**
	 * 
	 * @param A A[2][3][len]
	 * @return
	 */
	public static double[] determinant23(double[][][] A) {
		int len = A[0][0].length;
		double[] rlt = new double[len];
		for(int i=0;i<len;i++) {
			double f0 = A[0][0][i];
			double f1 = A[0][1][i];
			double f2 = A[0][2][i];
			double f3 = A[1][0][i];
			double f4 = A[1][1][i];
			double f5 = A[1][2][i];
			rlt[i] = Math.sqrt( (f1*f5-f4*f2)*(f1*f5-f4*f2) + (f0*f5-f3*f2)*(f0*f5-f3*f2) + (f0*f4-f3*f1)*(f0*f4-f3*f1) );
		}
		return rlt;
	}
	
	/**
	 * Golden Section Search for minimum value of f on interval [a,b] with accuracy eps
	 * 
	 * @param a
	 * @param b
	 * @param eps
	 * @param f
	 * @param debug
	 * @return
	 */
	public static double GoldenSectionSearch(double a, double b, double eps, MathFunc f, boolean debug) {
		if(debug) System.out.println("Golden Section Search "+f+" :");
		double r = (Math.sqrt(5)-1.0)/2.0; //0.618...
		double x = a + r*(b-a);
		double y = a + r*r*(b-a);
		double u = f.apply(new Variable(x));
		double v = f.apply(new Variable(y));
		while(b-a > eps) {
			if(debug) System.out.println(String.format("GSS Interval [%1.3f, %1.3f]: %s=%1.5f", a,b,f,u));
			if(u > v) {
				b = x;
				x = y;
				u = v;
				y = a + r*r*(b-a);
				v = f.apply(new Variable(y));
			} else {
				a = y;
				y = x;
				v = u;
				x = a + r*(b-a);
				u = f.apply(new Variable(x));
			}
		}
		return x;
	}
	
	public static void setDirichlet(Matrix stiff, Vector load, int matIndex, double value) {
		int row = matIndex;
		int col = matIndex;
		stiff.set(row, col, 1.0);
		load.set(row,value);
		for(int r=1;r<=stiff.getRowDim();r++) {
			if(r != row) {
				load.add(r,-stiff.get(r, col)*value);
				stiff.set(r, col, 0.0);
			}
		}
		for(int c=1;c<=stiff.getColDim();c++) {
			if(c != col) {
				stiff.set(row, c, 0.0);
			}
		}
	}

	public static void imposeDirichletCondition(Matrix stiff, Vector load, FiniteElement fe, Mesh mesh, MathFunc diri) {
		ElementList eList = mesh.getElementList();
		for(int i=1;i<=eList.size();i++) {
			NodeList nodes = eList.at(i).nodes;
			for(int j=1; j<=nodes.size(); j++) {
				Node n = nodes.at(j);
				if(n.getNodeType() == NodeType.Dirichlet) {
					Variable v = Variable.createFrom(diri, n, n.globalIndex); //bugfix 11/27/2013 Variable.createFrom(diri, n, 0);
					setDirichlet(stiff, load, fe.getGlobalIndex(mesh, eList.at(i), j),diri.apply(v));
				}
			}
		}
	}
	
	@Deprecated
	public static void imposeDirichletCondition(Matrix stiff, Vector load, Mesh mesh, MathFunc diri) {
		ElementList eList = mesh.getElementList();
		for(int i=1;i<=eList.size();i++) {
			Element e = eList.at(i);
			DOFList DOFs = e.getAllDOFList(DOFOrder.NEFV);
			for(int j=1;j<=DOFs.size();j++) {
				DOF dof = DOFs.at(j);
				GeoEntity ge = dof.getOwner();
				if(ge instanceof Node) {
					Node n = (Node)ge;
					if(n.getNodeType() == NodeType.Dirichlet) {
						Variable v = Variable.createFrom(diri, n, n.globalIndex); //bugfix 11/27/2013 Variable.createFrom(diri, n, 0);
						setDirichlet(stiff, load, dof.getGlobalIndex(),diri.apply(v));
					}
				} else if(ge instanceof EdgeLocal) {
					//2D单元（面）其中的局部边上的自由度
					EdgeLocal edge = (EdgeLocal)ge;
					if(edge.getBorderType() == NodeType.Dirichlet) {
						//TODO 以边的那个顶点取值？中点？
						//Variable v = Variable.createFrom(fdiri, ?, 0);
					}
					
				} else if(ge instanceof FaceLocal) {
					//3D单元（体）其中的局部面上的自由度
					FaceLocal face = (FaceLocal)ge;
					if(face.getBorderType() == NodeType.Dirichlet) {
						//TODO
					}
				} else if(ge instanceof Edge) {
					//1D单元（线段）上的自由度，其Dirichlet边界用结点来计算推出，而不需要专门标记单元
					VertexList vs = ((GeoEntity2D) ge).getVertices();
					for(int k=1;k<=vs.size();k++) {
						Node n = vs.at(k).globalNode();
						if(NodeType.Dirichlet == n.getNodeType()) {
							Variable v = Variable.createFrom(diri, n, 0);
							setDirichlet(stiff, load, dof.getGlobalIndex(),diri.apply(v));
						}
					}
				} else if(ge instanceof Face) {
					//2D单元（面）上的自由度，其Dirichlet边界用结点来计算推出，而不需要专门标记单元
					
					VertexList vs = ((GeoEntity2D) ge).getVertices();
					for(int k=1;k<=vs.size();k++) {
						Node n = vs.at(k).globalNode();
						if(NodeType.Dirichlet == n.getNodeType()) {
							Variable v = Variable.createFrom(diri, n, 0);
							setDirichlet(stiff, load, dof.getGlobalIndex(),diri.apply(v));
						}
					}
				} else if(ge instanceof Volume) {
					//3D单元（体）上的自由度，其Dirichlet边界用结点来计算推出，而不需要专门标记单元
					VertexList vs = ((GeoEntity3D) ge).getVertices();
					for(int k=1;k<=vs.size();k++) {
						Node n = vs.at(k).globalNode();
						if(NodeType.Dirichlet == n.getNodeType()) {
							Variable v = Variable.createFrom(diri, n, 0);
							setDirichlet(stiff, load, dof.getGlobalIndex(),diri.apply(v));
						}
					}
				}
			}
		}
	}

}
