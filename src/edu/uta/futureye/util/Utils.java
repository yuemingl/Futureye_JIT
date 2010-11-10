package edu.uta.futureye.util;

import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Map.Entry;

import edu.uta.futureye.algebra.Vector;
import edu.uta.futureye.core.DOF;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.core.NodeRefined;
import edu.uta.futureye.core.intf.Point;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.FAbstract;
import edu.uta.futureye.function.basic.FConstant;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.function.intf.FunctionDerivable;
import edu.uta.futureye.function.operator.FOBasic;
import edu.uta.futureye.function.shape.SFQuadraticLocal1D;

public class Utils {
	
	public static List<String> mergeList(List<String> a, List<String> b) {
		Set<String> set = new LinkedHashSet<String>();
		//TODO ??? 是否要做判断？
		if(a != null)
			set.addAll(a);
		if(b != null)
			set.addAll(b);
		List<String> rlt = new LinkedList<String>();
		rlt.addAll(set);
		return rlt;
	}
	
	public static double computeAngle2D(Point a1,Point a2,Point b1,Point b2) {
		Vector v1 = new Vector(2);
		v1.set(1, a2.coord(1)-a1.coord(1));
		v1.set(2, a2.coord(2)-a1.coord(2));
		
		Vector v2 = new Vector(2);
		v2.set(1, b2.coord(1)-b1.coord(1));
		v2.set(2, b2.coord(2)-b1.coord(2));
		
		//if(v1.norm2()<Constant.eps || v2.norm2()<Constant.eps), Math.acos(v) will be NaN

		double v = v1.dot(v2)/(v1.norm2()*v2.norm2());
		if(v > 1.0) v = 1.0;
		else if(v < -1.0) v = -1.0;
		return Math.acos(v);
		
	}

	public static double linearInterpolate(double x1, double x2, double x, 
			double y1, double y2) {
		double k = (y2 -y1)/(x2 -x1);
		double r = k*(x-x1)+y1;
		return r;
	}
	
	
	public static double linearInterpolate(Point p1, Point p2, Point p, 
			double y1, double y2) {
		double x1,x2,x;
		for(int i=1;i<=p1.dim();i++) {
			x1 = p1.coord(i);
			x2 = p2.coord(i);
			x = p.coord(i);
			if(Math.abs(x1-x2)>Constant.eps) {
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
	 * @param mesh
	 * @param u
	 * @param neighborBand
	 * @param weight
	 * @return
	 */
	public static Vector gaussSmooth(Mesh mesh, Vector u, int neighborBand, double weight) {
		NodeList list = mesh.getNodeList();
		int nNode = list.size();
		if(nNode != u.getDim()) {
			FutureEyeException e = new FutureEyeException("Node number of mesh != length of vector u: "+
					"nNode="+nNode+"  dim(u)="+u.getDim());
			e.printStackTrace();
			return null;
		}
		
	    Vector su = new Vector(nNode);
	    
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
		    	if(nbList.size() == 0) {
					FutureEyeException e = new FutureEyeException("No beighbors of Node "+node.globalIndex+", call mesh.computeNeiborNode() first!");
					e.printStackTrace();
					return null;
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
		    	if(nbList.size() == 0) {
					FutureEyeException e = new FutureEyeException("No beighbors of Node "+node.globalIndex+", call mesh.computeNeiborNode() first!");
					e.printStackTrace();
					return null;
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
	
	public static Function interplateFunctionOnElement(Function fun, Element e) {
		if(fun instanceof FConstant)
			return fun;
		Function rlt = new FConstant(0.0);
		int nNode = e.nodes.size();
		for(int i=1;i<=nNode;i++) {
			DOFList dofListI = e.getDOFList(i);
			for(int k=1;k<=dofListI.size();k++) {
				DOF dofI = dofListI.at(k);
				Variable var = Variable.createFrom(fun, dofI.getOwnerNode(), dofI.getGlobalNumber());
				FunctionDerivable PValue = new FConstant(fun.value(var));
				rlt = FOBasic.Plus(rlt, FOBasic.Mult(PValue, dofI.getShapeFunction()));
			}
		}
		return rlt;
	}
	
	//只是用于三角形线性元
	public static Map<String, FunctionDerivable> getFunctionComposeMap(Element e) {
		final Element fe =e;
		Map<String, FunctionDerivable> fInners = new HashMap<String, FunctionDerivable>();
		final List<String> varNamesInner = new LinkedList<String>();
		varNamesInner.add("r");
		varNamesInner.add("s");
		varNamesInner.add("t");
		fInners.put("x", new FAbstract(varNamesInner) {	
			@Override
			public double value(Variable v) {
				double rlt = 0.0;
				for(int i=1;i<=fe.nodes.size();i++)
					rlt += fe.nodes.at(i).coord(1)*v.get(varNamesInner.get(i-1));
				return rlt;
				
			}
		});
		fInners.put("y", new FAbstract(varNamesInner) {	
			@Override
			public double value(Variable v) {
				double rlt = 0.0;
				for(int i=1;i<=fe.nodes.size();i++)
					rlt += fe.nodes.at(i).coord(2)*v.get(varNamesInner.get(i-1));
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
		if(Math.abs(rlt-Math.PI) < Constant.eps) {
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
		if(Math.abs(rlt-Math.PI) < Constant.eps) {
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
		if(Math.abs(rlt-Math.PI) < Constant.eps ||
				Math.abs(rlt) < Constant.eps) {
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
	
}
