package edu.uta.futureye.util;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import edu.uta.futureye.algebra.SpaceVector;
import edu.uta.futureye.algebra.SparseVectorHashMap;
import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.SparseVector;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.algebra.solver.external.SolverJBLAS;
import edu.uta.futureye.core.DOF;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.core.NodeRefined;
import edu.uta.futureye.core.NodeType;
import edu.uta.futureye.core.geometry.Point;
import edu.uta.futureye.core.intf.AssemblerOld;
import edu.uta.futureye.function.FMath;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.FXY;
import edu.uta.futureye.function.basic.Vector2MathFunc;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.function.intf.ScalarShapeFunction;
import edu.uta.futureye.io.MatlabMatFileReader;
import edu.uta.futureye.io.MeshWriter;
import edu.uta.futureye.lib.assembler.AssemblerScalar;
import edu.uta.futureye.lib.shapefun.SFBilinearLocal2D;
import edu.uta.futureye.lib.shapefun.SFLinearLocal2D;
import edu.uta.futureye.lib.weakform.WeakFormDerivative;
import edu.uta.futureye.lib.weakform.WeakFormL22D;
import edu.uta.futureye.util.Constant;
import edu.uta.futureye.util.FutureyeException;
import edu.uta.futureye.util.PairDoubleInteger;
import edu.uta.futureye.util.Utils;
import edu.uta.futureye.util.container.ElementList;
import edu.uta.futureye.util.container.NodeList;

public class Tools {
	public static Vector computeDerivative(Mesh mesh, Vector U, String varName) {
		return computeDerivative(mesh, U, varName, 0.0);
		
		//return computeDerivativeFast(mesh, U, varName);
		
	}
	
	public static Vector computeDerivative(Mesh mesh, Vector U, String varName, double stableFactor) {
		//2011-06-29
		//没必要清除边界标记
		//mesh.clearBorderNodeMark();
		
		WeakFormDerivative weakForm = new WeakFormDerivative(varName);
		//stabilize
		weakForm.setParam(new Vector2MathFunc(U),stableFactor);
		
		AssemblerScalar assembler = new AssemblerScalar(mesh, weakForm);
		assembler.assemble();
		Matrix stiff = assembler.getStiffnessMatrix();
		Vector load = assembler.getLoadVector();
		
		//Solver solver = new Solver();
		//Vector w = solver.solveCGS(stiff, load);
		
		SolverJBLAS solver = new SolverJBLAS();
		Vector w = solver.solveDGESV(stiff, load);
		
		//
		//包含handing结点的网格计算导数的时候会产生间断，这里强制连续
		//
//		NodeList nodes = mesh.getNodeList();
//		for(int i=1;i<=nodes.size();i++) {
//			Node node = nodes.at(i);
//			if(node instanceof NodeRefined) {
//				NodeRefined nRefined = (NodeRefined)node;
//				if(nRefined.isHangingNode()) {
//					w.set(node.globalIndex,
//							w.get(nRefined.constrainNodes.at(1).globalIndex)*0.5+
//							w.get(nRefined.constrainNodes.at(2).globalIndex)*0.5);
//				}
//			}
//		}
		return w;
	}
	
	
	public static Vector computeDerivativeFast(Mesh mesh, Vector U, String varName) {
		Map<Integer, Double[]> map = new HashMap<Integer, Double[]>();
		ElementList eList = mesh.getElementList();
		Vector rltVector = U.copy();
		rltVector.setAll(0.0);
		
		MathFunc fU = new Vector2MathFunc(U);
		for(int i=1;i<=eList.size();i++) {
			Element e = eList.at(i);
			int N = e.nodes.size();
			double[] f = new double[N];
			for(int j=1;j<=N;j++) {
				Node node = e.nodes.at(j);
				Variable var = Variable.createFrom(fU, node, node.globalIndex);
				f[j-1] = fU.apply(var);
			}
			double[] a = Utils.computeBilinearFunctionCoef(
					e.nodes.toArray(new Point[0]), f);
			Double[] aa = new Double[N];
			for(int j=0;j<N;j++) 
				aa[j] = a[j];
			if(map.get(e.globalIndex) != null) {
				throw new FutureyeException("Duplicate element in mesh: i="+i+", e.globalIndex="+e.globalIndex);
			}
			map.put(e.globalIndex, aa);
		}
		
		NodeList nodes = mesh.getNodeList();
		for(int i=1;i<=nodes.size();i++) {
			Node node = nodes.at(i);
			//找到结点所属的单元，在每个单元上计算导数，最后取平均值
			ElementList list = node.belongToElements;
			if(list == null || list.size()==0) {
				throw new FutureyeException("node.belongToElements is empty!");
			}
			double dv = 0.0;
			int nE = list.size();
			
			//判断一个结点周围单元不包含Hanging Node的单元总数：
			//  如果大于零则利用这些单元计算导数
			//  如果等于零，则令包含1个Handing Node单元也参与导数计算
			//这样处理的结果会使导数更加光滑
			int threthHode = 0;
			int nTotal=0;
			for(int j=1;j<=nE;j++) {
				Element e = list.at(j);
				if(e.getHangingNode().size()==0) nTotal++;
			}
			if(nTotal==0) threthHode = 1;
			//threthHode = 4;//效果要好一些？
			nTotal = 0;
			for(int j=1;j<=nE;j++) {
				Element e = list.at(j);
				Double[] aa = map.get(e.globalIndex);
				MathFunc dd = null;
				//如果单元包含2个Handing Node，则不参与导数计算，这样处理的结果会使导数更加光滑
				if(e.getHangingNode().size()<=threthHode) {
					if(varName.equals("x")) {
						dd = new FXY(0.0,aa[3],aa[1]);//d(a1 + a2*x + a3*y + a4*x*y)/dx
					} else if(varName.equals("y")) {
						dd = new FXY(aa[3],0.0,aa[2]);//d(a1 + a2*x + a3*y + a4*x*y)/dy
					} else
						throw new FutureyeException("Parameter varName(="+varName+") should be x or y!");
					dv += dd.apply(new Variable("x",node.coord(1)).set("y",node.coord(2)));
					nTotal++;
				}
			}
			dv /= nTotal;
			rltVector.set(node.globalIndex,dv);
		}
		
		
		//
		//包含handing结点的网格计算导数的时候会产生间断，这里强制连续
		//
		for(int i=1;i<=nodes.size();i++) {
			Node node = nodes.at(i);
			if(node instanceof NodeRefined) {
				NodeRefined nRefined = (NodeRefined)node;
				if(nRefined.isHangingNode()) {
					rltVector.set(node.globalIndex,
							rltVector.get(nRefined.constrainNodes.at(1).globalIndex)*0.5+
							rltVector.get(nRefined.constrainNodes.at(2).globalIndex)*0.5);
				}
			}
		}
		return rltVector;
	}
	
	public static Vector computeLaplace2D(Mesh mesh, Vector U) {
		Vector ux = Tools.computeDerivative(mesh,U,"x");
		Vector uy = Tools.computeDerivative(mesh,U,"y");
		Vector uxx = Tools.computeDerivative(mesh,ux,"x");
		Vector uyy = Tools.computeDerivative(mesh,uy,"y");
		Vector LpU = FMath.axpy(1.0, uxx, uyy);
		return LpU;
	}
	
	public static void plotVector(Mesh mesh, String outputFolder, String fileName, Vector v, Vector ...vs) {
	    MeshWriter writer = new MeshWriter(mesh);
	    if(!outputFolder.isEmpty()) {
		    File file = new File("./"+outputFolder);
			if(!file.exists()) {
				file.mkdirs();
			}
	    }
	    writer.writeTechplot("./"+outputFolder+"/"+fileName, v, vs);
	}

	public static void plotFunction(Mesh mesh, String outputFolder, String fileName, MathFunc fun, MathFunc ...funs) {
//	    NodeList list = mesh.getNodeList();
//	    int nNode = list.size();
//		Variable var = new Variable();
//		Vector v = new SparseVector(nNode);
//	    for(int i=1;i<=nNode;i++) {
//	    	Node node = list.at(i);
//	    	var.setIndex(node.globalIndex);
//	    	for(int j=1;j<=node.dim();j++) {
//	    		if(fun.varNames().size()==node.dim())
//	    			var.set(fun.varNames().get(j-1), node.coord(j));
//	    	}
//	    	v.set(i, fun.value(var));
//	    }
		Vector v = Tools.function2vector(mesh, fun);
		if(funs.length > 0) {
			Vector[] vs = new SparseVector[funs.length];
			for(int i=0;i<funs.length;i++) {
				vs[i] = Tools.function2vector(mesh, funs[i]);
			}
		    plotVector(mesh,outputFolder,fileName,v,vs);
			
		} else {
			plotVector(mesh,outputFolder,fileName,v);
		}
	}
	
	public static Vector extendData(Mesh meshFrom, Mesh meshTo, Vector u, double deaultValue) {
		NodeList nodeTo = meshTo.getNodeList();
		int dimTo = nodeTo.size();
		Vector rlt = u.copy();
		rlt.setDim(dimTo);
		for(int i=1;i<=nodeTo.size();i++) {
			Node node = meshFrom.findNode(nodeTo.at(i));
			if(node != null) {
				rlt.set(i, u.get(node.globalIndex));
			} else {
				rlt.set(i, deaultValue);
			}
		}
		return rlt;
	}	
	/**
	 * 
	 * @param mesh1
	 * @param mesh2
	 * @param u
	 * @return
	 */
	public static Vector extractData(Mesh meshFrom, Mesh meshTo, Vector u) {
		NodeList nodeTo = meshTo.getNodeList();
		int dimTo = nodeTo.size();
		Vector rlt = u.copy();
		rlt.setDim(dimTo);
		for(int i=1;i<=nodeTo.size();i++) {
			Node node = meshFrom.findNode(nodeTo.at(i));
			if(node != null) {
				rlt.set(i, u.get(node.globalIndex));
			}
		}
		return rlt;
	}	
	
	/**
	 * Solve the L2 problem for a(x)
	 * 
	 * (k*\nabla{U},\phi) + (a(x)*U,\phi) = f
	 * 
	 * @param mesh
	 * @param U
	 * @param f
	 * @param k
	 * @return a(x)
	 */
	public static Vector solveParamInverse(Mesh mesh, Vector U, MathFunc f,MathFunc k,MathFunc diri) {
		HashMap<NodeType, MathFunc> mapNTF2 = new HashMap<NodeType, MathFunc>();
		mapNTF2.put(NodeType.Dirichlet, null);
		mesh.clearBorderNodeMark();
		mesh.markBorderNode(mapNTF2);
		
		//Weak form
		WeakFormL22D weakFormL2 = new WeakFormL22D();
		//Right hand side
		weakFormL2.setF(f);
		
		
		//Parameters
//利用弱形式内部定义的单元上的分片导数计算Ux,Uy		
		weakFormL2.setParam(
				k, new Vector2MathFunc(U)
			);

//求解方程计算导数
//		Vector Ux = Tools.computeDerivative(mesh, U, "x",0.0);
//		Vector Uy = Tools.computeDerivative(mesh, U, "y",0.0);
		
//快速计算导数		
//		Vector Ux = Tools.computeDerivativeFast(mesh, U, "x");
//		Vector Uy = Tools.computeDerivativeFast(mesh, U, "y");
//		weakFormL2.setParam(
//				k, new Vector2Function(U), 
//				new Vector2Function(Ux),
//				new Vector2Function(Uy)
//			);		
		
		AssemblerScalar assembler = new AssemblerScalar(mesh, weakFormL2);
		//System.out.println("Begin Assemble...solveParamInverse");
		assembler.printInfo(false);
		assembler.assemble();
		Matrix stiff = assembler.getStiffnessMatrix();
		Vector load = assembler.getLoadVector();
		assembler.imposeDirichletCondition(diri);
		//System.out.println("Assemble done!");
		
		//Tools.plotVector(mesh, "", "load2.dat", load);
		
		SolverJBLAS solver = new SolverJBLAS();
		Vector u = solver.solveDGESV(stiff, load);
		return u;
	}
	
	public static Vector function2vector(Mesh mesh, MathFunc f) {
		NodeList nodes = mesh.getNodeList();
		Vector v = new SparseVectorHashMap(nodes.size());
		for(int j=1;j<=nodes.size();j++) {
			v.set(j, f.apply(Variable.createFrom(f, nodes.at(j), j)));
		}
		return v;
	}
	
	public static ElementList computeRefineElement(Mesh mesh, Vector v, double persent) {
	    Vector v_smooth = Utils.gaussSmooth(mesh, v, 1, 0.5);

	    mesh.computeNodeBelongsToElements();
	    mesh.computeGlobalEdge();
	    mesh.computeNeighborElements();
	    
		ElementList eList = mesh.getElementList();
		ElementList eToRefine = new ElementList();

		List<PairDoubleInteger> list = new ArrayList<PairDoubleInteger>();
		for(int i=1;i<=eList.size();i++) {
			Element e = eList.at(i);
			//L2Norm(v-v_smooth) on element e
			Vector v1 = new SparseVectorHashMap(e.nodes.size());
			for(int j=1;j<=e.nodes.size();j++) {
				Node node = e.nodes.at(j);
				v1.set(j, v.get(node.globalIndex)-v_smooth.get(node.globalIndex));
			}
			list.add(new PairDoubleInteger(v1.norm2(),i));
		}
		Collections.sort(list, new Comparator<PairDoubleInteger>() {
			@Override
			public int compare(PairDoubleInteger o1, PairDoubleInteger o2) {
				return o1.d < o2.d ? 1 : -1;
			}
		});
		int nThreshold = (int) Math.floor(persent*eList.size());
		for(int i=1;i<=nThreshold;i++) {
			eToRefine.add(eList.at(list.get(i).i));
		}
		
		//如果一个单元相邻单元都细化了，该单元自动细化
		for(int i=1;i<=eList.size();i++) {
			Element e = eList.at(i);
			ElementList eNeighbors = e.neighbors;
			int findTimes = 0;
			for(int k=1;k<=eNeighbors.size();k++) {
				Element en = eNeighbors.at(k);
				for(int j=1;j<=eToRefine.size();j++) {
					if(eToRefine.at(j).globalIndex == en.globalIndex) {
						findTimes++;
						break;
					}
				}
			}
			if(findTimes == eNeighbors.size()) {
				boolean find = false;
				for(int j=1;j<=eToRefine.size();j++) {
					if(eToRefine.at(j).globalIndex == e.globalIndex) {
						find = true;
						break;
					}
				}
				if(!find)
					eToRefine.add(e);
			}
		}
		
		return eToRefine;
	}	
	
	public static ElementList computeRefineElementMax(Mesh mesh, 
			Vector indicator, double persent) {

	    mesh.computeNodeBelongsToElements();
	    mesh.computeGlobalEdge();
	    mesh.computeNeighborElements();
	    
		ElementList eList = mesh.getElementList();
		ElementList eToRefine = new ElementList();

		List<PairDoubleInteger> list = new ArrayList<PairDoubleInteger>();
		for(int i=1;i<=eList.size();i++) {
			Element e = eList.at(i);
			//sum(indicator) on element e
			Vector v1 = new SparseVectorHashMap(e.nodes.size());
			for(int j=1;j<=e.nodes.size();j++) {
				v1.set(j, indicator.get(e.nodes.at(j).globalIndex));
			}
			list.add(new PairDoubleInteger(FMath.sum(v1),i));
		}
		Collections.sort(list, new Comparator<PairDoubleInteger>() {
			@Override
			public int compare(PairDoubleInteger o1, PairDoubleInteger o2) {
				return o1.d < o2.d ? 1 : -1;
			}
		});
		int nThreshold = (int) Math.floor(persent*eList.size());
		for(int i=1;i<=nThreshold;i++) {
			eToRefine.add(eList.at(list.get(i).i));
		}
		
		//如果一个单元相邻单元都细化了，该单元自动细化
		for(int i=1;i<=eList.size();i++) {
			Element e = eList.at(i);
			ElementList eNeighbors = e.neighbors;
			int findTimes = 0;
			for(int k=1;k<=eNeighbors.size();k++) {
				Element en = eNeighbors.at(k);
				for(int j=1;j<=eToRefine.size();j++) {
					if(eToRefine.at(j).globalIndex == en.globalIndex) {
						findTimes++;
						break;
					}
				}
			}
			if(findTimes == eNeighbors.size()) {
				boolean find = false;
				for(int j=1;j<=eToRefine.size();j++) {
					if(eToRefine.at(j).globalIndex == e.globalIndex) {
						find = true;
						break;
					}
				}
				if(!find)
					eToRefine.add(e);
			}
		}
		
		return eToRefine;
	}	
	
	public static void assignLinearShapFunction(Mesh mesh) {
		ScalarShapeFunction[] shapeFun = null;
		ScalarShapeFunction[] shapeFunHalf = null;
		ScalarShapeFunction[] shapeFunRect = null;
		ScalarShapeFunction[] shapeFunRectHalf = null;
		
		//Assign degree of freedom to element
		shapeFun = new SFLinearLocal2D[3];
		shapeFunHalf = new SFLinearLocal2D[3];
		for(int i=0;i<3;i++) {
			shapeFun[i] = new SFLinearLocal2D(i+1);
			shapeFunHalf[i] = new SFLinearLocal2D(i+1,0.5);
		}
		
		shapeFunRect = new SFBilinearLocal2D[4];
		shapeFunRectHalf = new SFBilinearLocal2D[4];
		for(int i=0;i<4;i++) {
			shapeFunRect[i] = new SFBilinearLocal2D(i+1);
			shapeFunRectHalf[i] = new SFBilinearLocal2D(i+1,0.5);
		}
		
		//Assign shape function to DOF
		for(int i=1;i<=mesh.getElementList().size();i++) {
			Element e = mesh.getElementList().at(i);
			e.clearAllDOF();
			if(e.nodes.size() == 4) {
				//Asign degree of freedom to element
				int nDofLocalIndexCounter = 0;
				for(int j=1;j<=e.nodes.size();j++) {
					//Asign shape function to DOF
					if(e.nodes.at(j) instanceof NodeRefined) {
						NodeRefined nRefined = (NodeRefined)e.nodes.at(j);
						if(nRefined.isHangingNode()) {
							DOF dof  = new DOF(++nDofLocalIndexCounter,nRefined.constrainNodes.at(1).globalIndex,
									shapeFunRectHalf[j-1]);
							e.addNodeDOF(j, dof);
							DOF dof2 = new DOF(++nDofLocalIndexCounter,nRefined.constrainNodes.at(2).globalIndex,
									shapeFunRectHalf[j-1]);
							e.addNodeDOF(j, dof2);
						} else {
							DOF dof = new DOF(++nDofLocalIndexCounter,e.nodes.at(j).globalIndex,shapeFunRect[j-1]);
							e.addNodeDOF(j, dof);				
						}
					} else {
						DOF dof = new DOF(++nDofLocalIndexCounter,e.nodes.at(j).globalIndex,shapeFunRect[j-1]);
						e.addNodeDOF(j, dof);
					}
				}
			} else if(e.nodes.size() == 3) {
				int nDofLocalIndexCounter = 0;
				for(int j=1;j<=e.nodes.size();j++) {
					//Asign shape function to DOF
					if(e.nodes.at(j) instanceof NodeRefined) {
						NodeRefined nRefined = (NodeRefined)e.nodes.at(j);
						if(nRefined.isHangingNode()) {
							DOF dof  = new DOF(++nDofLocalIndexCounter,nRefined.constrainNodes.at(1).globalIndex,
									shapeFunHalf[j-1]);
							e.addNodeDOF(j, dof);
							DOF dof2 = new DOF(++nDofLocalIndexCounter,nRefined.constrainNodes.at(2).globalIndex,
									shapeFunHalf[j-1]);
							e.addNodeDOF(j, dof2);
						} else {
							DOF dof = new DOF(++nDofLocalIndexCounter,e.nodes.at(j).globalIndex,shapeFun[j-1]);
							e.addNodeDOF(j, dof);				
						}
					} else {
						DOF dof = new DOF(++nDofLocalIndexCounter,e.nodes.at(j).globalIndex,shapeFun[j-1]);
						e.addNodeDOF(j, dof);
					}
				}
				
			} else {
				System.out.println("Error: e.nodes.size()="+e.nodes.size());
			}
		}
	}

	public static Vector interplateFrom(Mesh oldMesh, Mesh newMesh, Vector2MathFunc vecFun) {
		int nNode = newMesh.getNodeList().size();
		NodeList nodes = newMesh.getNodeList();
		Vector rlt = new SparseVectorHashMap(nNode);
		for(int i=1;i<=nNode;i++) {
			Variable v = new Variable();
			v.set("x",nodes.at(i).coord(1));
			v.set("y",nodes.at(i).coord(2));
			double val = vecFun.apply(v);
			rlt.set(i,val);
		}
		return rlt;
	}
	
	public static void constrainHangingNodes(Mesh mesh, Vector v) {
		NodeList nodes = mesh.getNodeList();
		for(int i=1;i<=nodes.size();i++) {
			Node node = nodes.at(i);
			if(node instanceof NodeRefined) {
				NodeRefined nRefined = (NodeRefined)node;
				if(nRefined.isHangingNode()) {
					v.set(node.globalIndex,
							v.get(nRefined.constrainNodes.at(1).globalIndex)*0.5+
							v.get(nRefined.constrainNodes.at(2).globalIndex)*0.5);
				}
			}
		}
	}
	
	/**
	 * Read values of f(x,y), where (x,y) defined on <tt>mesh</tt>,
	 * return a vector defined on <tt>mesh</tt> containing the values. 
	 * 
	 * @param matFile
	 * @param mesh
	 * @param xName x: 2D matlab matrix
	 * @param yName y: 2D matlab matrix
	 * @param vName f(x,y): 2D matlab matrix 
	 * @return
	 */
	public static Vector read2DFunctionValues(Mesh mesh, String matFile, String xName, String yName, String vName) {
		MatlabMatFileReader r = new MatlabMatFileReader(matFile);
		Matrix mx = r.getMatrix(xName);
		Matrix my = r.getMatrix(yName);
		Matrix mv = r.getMatrix(vName);
		Vector rlt = new SpaceVector(mesh.getNodeList().size());
		int m = mv.getRowDim();
		int n = mv.getColDim();
		for(int i=1;i<=m;i++) {
			for(int j=1;j<=n;j++) {
				double x = mx.get(i, j);
				double y = my.get(i, j);
				double v = mv.get(i, j);
				Node node = mesh.findNode(new Node(0,x,y));
				if(node != null) {
					rlt.set(node.globalIndex, v);
				}
			}
		}
		return rlt;
	}
	
	/**
	 * 
	 * @param mesh
	 * @param x
	 * @param sortAsc
	 * @return
	 */
	public static NodeList getNodeListOnX(Mesh mesh, double x,final boolean sortAsc) {
		NodeList nodes = mesh.getNodeList();
		NodeList rlt = new NodeList();
		for(int i=1;i<=nodes.size();i++) {
			Node node = nodes.at(i);
			if(Math.abs(node.coord(1)-x)<Constant.meshEps) {
				rlt.add(node);
			}
		}
		List<Node> l = rlt.toList();
		Collections.sort(l, new Comparator<Node>() {
			@Override
			public int compare(Node o1, Node o2) {
				if(o1.coord(2)>o2.coord(2))
					return sortAsc?1:-1;
				else
					return sortAsc?-1:1;
			}
		});
		return rlt;
	}
	
	/**
	 * 
	 * @param mesh
	 * @param y
	 * @param sortAsc: sort by x
	 * @return
	 */
	public static NodeList getNodeListOnY(Mesh mesh, double y,final boolean sortAsc) {
		NodeList nodes = mesh.getNodeList();
		NodeList rlt = new NodeList();
		for(int i=1;i<=nodes.size();i++) {
			Node node = nodes.at(i);
			if(Math.abs(node.coord(2)-y)<Constant.meshEps) {
				rlt.add(node);
			}
		}
		List<Node> l = rlt.toList();
		Collections.sort(l, new Comparator<Node>() {
			@Override
			public int compare(Node o1, Node o2) {
				if(o1.coord(1)>o2.coord(1))
					return sortAsc?1:-1;
				else
					return sortAsc?-1:1;
			}
		});
		return rlt;
	}
	
	/**
	 * Scale the vectors <code>vs</code> according the maximum and minimum values
	 * of the first parameter <code>first</code> 
	 * @param first
	 * @param vs
	 */
	public static void scaleVectorsByFirst(Vector first, Vector... vs) {
		double max = FMath.max(first);
		double min = FMath.min(first);
		for(int i=0; i<vs.length; ++i) {
			double mx = FMath.max(vs[i]);
			double mn = FMath.min(vs[i]);
			vs[i].shift(-mn);
			double factor = (max-min)/(mx-mn);
			vs[i].scale(factor);
			vs[i].shift(mn);
		}
	}
	
	
}