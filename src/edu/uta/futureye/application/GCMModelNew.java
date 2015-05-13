package edu.uta.futureye.application;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;

import edu.uta.futureye.algebra.SparseVectorHashMap;
import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.algebra.solver.external.SolverJBLAS;
import edu.uta.futureye.application.MouseHead.TailType;
import edu.uta.futureye.core.DOF;
import edu.uta.futureye.core.EdgeLocal;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.core.NodeLocal;
import edu.uta.futureye.core.NodeRefined;
import edu.uta.futureye.core.NodeType;
import edu.uta.futureye.core.Refiner;
import edu.uta.futureye.core.Vertex;
import edu.uta.futureye.core.intf.Assembler;
import edu.uta.futureye.function.AbstractFunction;
import edu.uta.futureye.function.FMath;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.basic.FDelta;
import edu.uta.futureye.function.basic.Vector2Function;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.function.intf.ScalarShapeFunction;
import edu.uta.futureye.io.MeshReader;
import edu.uta.futureye.io.MeshWriter;
import edu.uta.futureye.lib.assembler.AssemblerScalar;
import edu.uta.futureye.lib.shapefun.SFBilinearLocal2D;
import edu.uta.futureye.lib.shapefun.SFLinearLocal2D;
import edu.uta.futureye.lib.shapefun.SFQuadraticLocal2DFast;
import edu.uta.futureye.lib.weakform.WeakFormL22D;
import edu.uta.futureye.lib.weakform.WeakFormLaplace2D;
import edu.uta.futureye.util.Constant;
import edu.uta.futureye.util.FutureyeException;
import edu.uta.futureye.util.PairDoubleInteger;
import edu.uta.futureye.util.Utils;
import edu.uta.futureye.util.container.ElementList;
import edu.uta.futureye.util.container.NodeList;
import edu.uta.futureye.util.container.ObjList;

/**
 *   
 * @author liuyueming
 */
public class GCMModelNew {
	String outputFolder = null;
	boolean debug = false;
	
	public GCMModelNew(String outputFolder) {
		this.outputFolder = outputFolder;
	}
	
	public void plotFunction(Mesh mesh, Function fun, String fileName) {
	    NodeList list = mesh.getNodeList();
	    int nNode = list.size();
		Variable var = new Variable();
		Vector v = new SparseVectorHashMap(nNode);
	    for(int i=1;i<=nNode;i++) {
	    	Node node = list.at(i);
	    	var.setIndex(node.globalIndex);
	    	var.set("x", node.coord(1));
	    	var.set("y", node.coord(2));
	    	v.set(i, fun.value(var));
	    }
	    plotVector(mesh,v,fileName);
	}
	
	public Vector discreteFunction(Mesh mesh,Function fun) {
	    NodeList list = mesh.getNodeList();
	    int nNode = list.size();
		Variable var = new Variable();
		Vector v = new SparseVectorHashMap(nNode);
	    for(int i=1;i<=nNode;i++) {
	    	Node node = list.at(i);
	    	var.setIndex(node.globalIndex);
	    	var.set("x", node.coord(1));
	    	var.set("y", node.coord(2));
	    	v.set(i, fun.value(var));
	    }
	    return v;
	}
	
	public void plotVector(Mesh mesh, Vector v, String fileName) {
	    MeshWriter writer = new MeshWriter(mesh);
	    if(!outputFolder.isEmpty()) {
		    File file = new File("./"+outputFolder);
			if(!file.exists()) {
				file.mkdirs();
			}
	    }
	    writer.writeTechplot("./"+outputFolder+"/"+fileName, v);
	}
	
	Vector getGCMCoef(double sn_1, double sn) {
		Vector A = new SparseVectorHashMap(4);
		double I0 =  3*(sn_1-sn) - 2*sn_1*(Math.log(sn_1)-Math.log(sn));
		double A1 = (1.0/I0)*(
				(-3.0/2.0)*(Math.pow(sn_1, 4)-Math.pow(sn, 4))+
				(10.0/3.0)*sn_1*(Math.pow(sn_1, 3)-Math.pow(sn, 3))+
				(-2.0)*sn_1*sn_1*(Math.pow(sn_1, 2)-Math.pow(sn, 2))
				);
		double A2 = (1.0/I0)*(
				(-10.0/3.0)*(Math.pow(sn_1, 3)-Math.pow(sn, 3))+
				4.0*sn_1*(Math.pow(sn_1, 2)-Math.pow(sn, 2))
				);
		double A3 = (-2.0/I0)*(Math.log(sn_1)-Math.log(sn));
		double A4 = (-2.0/I0)*(Math.pow(sn_1, 2)-Math.pow(sn, 2));
		A.set(1, A1);
		A.set(2, A2);
		A.set(3, A3);
		A.set(4, A4);
		return A;
	}

	/**
	 * 
	 * @param mesh
	 * @param N
	 * @param s
	 * @param phi
	 * @param tailT 
	 */
	public  void solveGCM(Mesh mesh, int N, double[]s, Vector[]phi, Vector tailT) {
		int dim = tailT.getDim();
		Vector sumLaplaceQ = new SparseVectorHashMap(dim);
		Vector sumQ_x= new SparseVectorHashMap(dim);
		Vector sumQ_y= new SparseVectorHashMap(dim);
//		Vector sumLaplaceQ_real = new SparseVectorHashMap(dim);
//		Vector sumQ_x_real = new SparseVectorHashMap(dim);
//		Vector sumQ_y_real = new SparseVectorHashMap(dim);
		Vector[] q = new Vector[N];
		Vector T_laplace = Tools.computeLaplace2D(mesh, tailT);
		Vector T_x = Tools.computeDerivative(mesh, tailT, "x");
		Vector T_y = Tools.computeDerivative(mesh, tailT, "y");
		if(debug) {
			plotVector(mesh, T_laplace, "T_laplace.dat");
			plotVector(mesh, T_x, "T_x.dat");
			plotVector(mesh, T_y, "T_y.dat");
		}
		Vector q_0 = new SparseVectorHashMap(dim);
		for(int i=0;i<N-1;i++) {
			q[i] = solveGCM(mesh,i,
					s[i],s[i+1],
					sumLaplaceQ,sumQ_x,sumQ_y,
					//sumLaplaceQ_real,sumQ_x_real,sumQ_y_real,
					T_laplace,T_x,T_y,
					new Vector2Function(phi[i]),
					q_0,phi[i],
					false //flag, if linearize equation q_n
					);
		
			q_0 = q[i];
			
			//sum_{j=1}^{N-1} \Laplace{q}
			sumLaplaceQ.axpy(1.0, Tools.computeLaplace2D(mesh, q[i]));
			sumQ_x.axpy(1.0, Tools.computeDerivative(mesh, q[i], "x"));
			sumQ_y.axpy(1.0, Tools.computeDerivative(mesh, q[i], "y"));//bugfix add 2011-5-7
			if(debug) {
				plotVector(mesh, sumLaplaceQ, "GCM_sumLaplaceQ_"+i+".dat");
				plotVector(mesh, sumQ_x, "GCM_sumQ_x_"+i+".dat");
				plotVector(mesh, sumQ_y, "GCM_sumQ_y_"+i+".dat");
			}
			
			//!!!TEST!!! Real q[i] = phi[i]
//			sumLaplaceQ_real.axpy(1.0,  Tools.computeLaplace2D(mesh, phi[i]));
//			sumQ_x_real.axpy(1.0,  Tools.computeDerivative(mesh, phi[i], "x"));
//			sumQ_y_real.axpy(1.0,  Tools.computeDerivative(mesh, phi[i], "y"));
//			plotVector(mesh, sumLaplaceQ_real, "GCM_sumLaplaceQ_real_"+i+".dat");
//			plotVector(mesh, sumQ_x_real, "GCM_sumQ_x_real_"+i+".dat");
//			plotVector(mesh, sumQ_y_real, "GCM_sumQ_y_real_"+i+".dat");
			
		}
		
		//用最近处光源重构(注意this.delta的位置)：v_tidle -> u ->　a
		Vector v_tidle = new SparseVectorHashMap(dim);
		for(int i=0;i<N-1;i++) {
			double h = s[i] - s[i+1];
			v_tidle.add(-s[N-1]*s[N-1]*h, q[i]);
			plotVector(mesh, v_tidle, String.format("v_tidle_N_%02d.dat",i));
		}
		v_tidle.add(s[N-1]*s[N-1], tailT);
		plotVector(mesh, v_tidle, "v_tidle_N.dat");
		
		Vector u = new SparseVectorHashMap(dim);
		for(int i=1;i<=dim;i++) {
			u.set(i, Math.pow(Math.E, v_tidle.get(i)));
		}
		plotVector(mesh, u, "u_N.dat");
		Math.ce
		
//		Vector a = solveParamInverse(mesh, u);
//		plotVector(mesh, a, "alpha_recon_N.dat");
		
		
//		//模拟计算出来的phi除了提供边界条件，区域内部的值也有，可以用来重构真实的a(x)
//		//*** v_tidle_real:= ln(u) ***
//		Vector v_tidle_real = new SparseVectorHashMap(dim);
//		for(int i=0;i<N-1;i++) {
//			double h = s[i] - s[i+1];
//			v_tidle_real.add(-s[N-1]*s[N-1]*h, phi[i]);
//			plotVector(mesh, v_tidle_real, String.format("v_tidle_N_real_%02d.dat",i));
//		}
//		v_tidle_real.add(s[N-1]*s[N-1], tailT);
//		plotVector(mesh, v_tidle_real, "v_tidle_N_real.dat");
//		
//		Vector u_real = new SparseVectorHashMap(dim);
//		for(int i=1;i<=dim;i++) {
//			u_real.set(i, Math.pow(Math.E, v_tidle_real.get(i)));
//		}
//		plotVector(mesh, u_real, "u_N_real.dat");
//		
//		Vector a_real = solveParamInverse(mesh, u_real);
//		plotVector(mesh, a_real, "alpha_recon_N_real.dat");
		
	}
	
	public Vector solveGCM(Mesh mesh,int n_1,
			double sn_1, double sn,
			Vector sumLaplaceQ,
			Vector sumQ_x, Vector sumQ_y,
			Vector laplaceT,
			Vector T_x, Vector T_y,
			Function diri,
			Vector q_0,
			Vector phi, //for check only
			boolean linearize
			) {
		
		//Mark border type
		HashMap<NodeType, Function> mapNTF = new HashMap<NodeType, Function>();
		mapNTF.clear();
		mapNTF.put(NodeType.Dirichlet, null);
//		mapNTF.put(NodeType.Dirichlet, new FAbstract("x","y"){
//			@Override
//			public double value(Variable v) {
//				double y = v.get("y");
//				if(Math.abs(y - 3.0) < Constant.eps)
//					return 1.0;
//				else
//					return -1.0;
//			}
//		});
//		mapNTF.put(NodeType.Robin, null);
		
		double h = sn_1 - sn; //错误的写成：sn-1 - sn，导致结果异常，以后变量名称一定不要这样写了
		Vector A = getGCMCoef(sn_1,sn);
		if(debug) {
			System.out.println("----------");
			for(int i=1;i<=A.getDim();i++) {
				System.out.print("A"+i+"="+A.get(i)+"  ");
			}
			System.out.println("\n----------");
		}
		//Right hand side of equation q_n
		int dim = sumLaplaceQ.getDim();
		Vector f1 = new SparseVectorHashMap(dim); //zero vector
		//??? + - ??? A.get(3)*h  //2011-5-7 "-" OK
		f1.add(-A.get(3)*h, sumLaplaceQ);
		f1.add(A.get(3), laplaceT);
		Vector f2x = new SparseVectorHashMap(dim); //zero vector
		f2x.add(h, sumQ_x);
		f2x.add(-1.0, T_x);
		f2x = FMath.axMuly(A.get(4), f2x, f2x);
		if(debug) plotVector(mesh, f2x, "GCM_f2x_"+n_1+".dat");
		Vector f2y = new SparseVectorHashMap(dim); //zero vector
		f2y.add(h, sumQ_y);
		f2y.add(-1.0, T_y);
		f2y = FMath.axMuly(A.get(4), f2y, f2y);
		if(debug) plotVector(mesh, f2y, "GCM_f2y_"+n_1+".dat");
		Vector f = FMath.axpy(1.0, f2x, f2y);
		if(debug) plotVector(mesh, f, "GCM_f2xpf2y_"+n_1+".dat");
		f.add(f1);
		if(debug) plotVector(mesh, f, "GCM_"+n_1+"_f.dat");
		
		//b1, b2
		Vector q_x_k_1 = null;
		Vector q_y_k_1 = null;
		
		//------------------------Nonlinear iteration-----------------------
		//Vector qn_1 = q_0;
		Vector qn_1 = new SparseVectorHashMap(dim);
		Vector qn = null;
		int nMaxIter = 15;
		if(linearize) nMaxIter = 1;
		for(int iter=1;iter<=nMaxIter;iter++) {
			
			//\nabla{q_n} term coefficient (include linearized nonlinear-term) 
			q_x_k_1 = Tools.computeDerivative(mesh, qn_1, "x");
			q_y_k_1 = Tools.computeDerivative(mesh, qn_1, "y");
			if(debug) {
				plotVector(mesh, q_x_k_1, "GCM_q_x_km1_"+n_1+"_"+iter+".dat");
				plotVector(mesh, q_y_k_1, "GCM_q_y_km1_"+n_1+"_"+iter+".dat");
			}
			Vector b1 = new SparseVectorHashMap(dim); //zero vector
			b1.add(A.get(2)*h, sumQ_x);
			b1.add(-A.get(2), T_x);
			if(!linearize) b1.add(-A.get(1), q_x_k_1);
			Vector b2 = new SparseVectorHashMap(dim); //zero vector
			b2.add(A.get(2)*h, sumQ_y);
			b2.add(-A.get(2), T_y);
			if(!linearize) b2.add(-A.get(1), q_y_k_1);
			if(debug) {
				plotVector(mesh, b1, "GCM_"+n_1+"_b1_"+iter+".dat");
				plotVector(mesh, b2, "GCM_"+n_1+"_b2_"+iter+".dat");
			}
			//check lhs with rhs
			if(debug && phi !=  null) {
				Vector phi_laplace = Tools.computeLaplace2D(mesh, phi);
				Vector phi_x = Tools.computeDerivative(mesh, phi, "x");
				Vector phi_y = Tools.computeDerivative(mesh, phi, "y");
				plotVector(mesh, phi_laplace, "GCM_"+n_1+"_phi_laplace_"+iter+".dat");
				plotVector(mesh, phi_x, "GCM_"+n_1+"_phi_x_"+iter+".dat");
				plotVector(mesh, phi_y, "GCM_"+n_1+"_phi_y_"+iter+".dat");
				Vector b1r = new SparseVectorHashMap(dim); //zero vector
				b1r.add(A.get(2)*h, sumQ_x);
				b1r.add(-A.get(2),  T_x);
				b1r.add(-A.get(1),  phi_x);
				Vector b2r = new SparseVectorHashMap(dim); //zero vector
				b2r.add(A.get(2)*h, sumQ_y);
				b2r.add(-A.get(2),  T_y);
				b2r.add(-A.get(1),  phi_y);
				Vector lhs_x = FMath.axMuly(1.0, phi_x, b1r);
				Vector lhs_y = FMath.axMuly(1.0, phi_y, b2r);
				Vector lhs = FMath.axpy(1.0, lhs_x, lhs_y);
				plotVector(mesh, lhs, "GCM_"+n_1+"_lhsx__lhsy_"+iter+".dat");
				lhs.add(phi_laplace);
				plotVector(mesh, lhs, "GCM_"+n_1+"_lhs_"+iter+".dat");
				//!!!TEST!!! Real \nabla{q_n} term coefficient 
				//b1 = b1r;
				//b2 = b2r;
			}
			
			WeakFormGCM weakForm = new WeakFormGCM();
			weakForm.setF(new Vector2Function(f));
			weakForm.setParam(
					new FC(-1.0),//-1 注意，有负号!!!
					new FC(0.0), //-\eps{q_n}!!! \eps-0.1
					new Vector2Function(b1),
					new Vector2Function(b2)
				);
			
			mesh.clearBorderNodeMark();
			mesh.markBorderNode(mapNTF);

			AssemblerScalar assembler = new AssemblerScalar(mesh, weakForm);
			System.out.println("Begin Assemble...solveGCM");
			assembler.assemble();
			Matrix stiff = assembler.getStiffnessMatrix();
			Vector load = assembler.getLoadVector();
			//Dirichlet condition
			assembler.imposeDirichletCondition(diri);
			System.out.println("Assemble done!");

			SolverJBLAS solver = new SolverJBLAS();
			qn = solver.solveDGESV(stiff, load);
			plotVector(mesh, qn, String.format("GCM_qn_%02d_%02d.dat",n_1,iter));
			
			if(debug) {
				Vector q_x_k = Tools.computeDerivative(mesh, qn, "x");
				Vector q_y_k = Tools.computeDerivative(mesh, qn, "y");
				Vector lhs_x = FMath.axMuly(1.0, q_x_k, b1);
				Vector lhs_y = FMath.axMuly(1.0, q_y_k, b2);
				Vector lhs = FMath.axpy(1.0, lhs_x, lhs_y);
				plotVector(mesh, lhs, "GCM_"+n_1+"_lhsx__lhsy_cmp_"+iter+".dat");
				Vector qn_laplace = Tools.computeLaplace2D(mesh, qn);
				lhs.add(qn_laplace);
				plotVector(mesh, lhs, "GCM_"+n_1+"_lhs_cmp_"+iter+".dat");
				Vector qn_equation = FMath.axpy(-1.0, f, lhs);
				plotVector(mesh, qn_equation, "GCM_"+n_1+"_qn_equation_cmp_"+iter+".dat");
			}
			
			//Stop criterion
			double eps = 0.01;
			System.out.println("GCM Iteration "+iter);
			if(iter > 1) {
				double norm = FMath.axpy(-1.0, qn_1, qn).norm2();
				if(norm < eps) {
					System.out.println("GCM norm < 1.0, STOP");
					break;
				}
			}
			qn_1 = qn;
		}
		return qn;
	}
	
	public static void sovleMouseHead() {
		GCMModel gcm = new GCMModel();
		gcm.outputFolder = "MouseHeadGCM";
		gcm.setDelta(0.0, 0.9);  //S13 y(0.8->0.9)
		//gcm.setDelta(1.2, 0.4);  //S9 x(1.1->1.2)

		//Read a triangle mesh from an input file
		//包围老鼠头的矩形规则区域
		//String gridName = "mouse";
		//String gridName = "mouse2";
		//String gridName = "mouse3";
		String gridName = "mouse4";
		
		MeshReader reader2 = new MeshReader(gridName+"_omega2.grd");
		Mesh meshOmega2 = reader2.read2DMesh();
		meshOmega2.computeNodeBelongsToElements();
		meshOmega2.computeNeighborNodes();

		int N = 3;
		Vector u1 = DataReader.readVector("MouseHead\\left_u1.dat");
		Vector u2 = DataReader.readVector("MouseHead\\left_u1.dat");
		Vector u3 = DataReader.readVector("MouseHead\\left_u1.dat");
		Vector[] phi = new Vector[3];
		phi[0] = u1;
		phi[1] = u2;
		phi[2] = u3;		
		double[] s = {1.0,1.1,1.2};

		Vector tailT = DataReader.readVector("MouseHead\\left_u1.dat");
		gcm.solveGCM(meshOmega2, N, s, phi, tailT);
		u1.print();		
	}
	
	public static void main(String[] args) {
		 sovleMouseHead();
	}
	
}
