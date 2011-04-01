package edu.uta.futureye.application;

import java.io.File;
import java.io.FileOutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.HashMap;

import edu.uta.futureye.algebra.Solver;
import edu.uta.futureye.algebra.SparseVector;
import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.core.NodeType;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.DuDn;
import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.basic.Vector2Function;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.function.operator.FMath;
import edu.uta.futureye.io.MeshReader;
import edu.uta.futureye.io.MeshWriter;
import edu.uta.futureye.lib.assembler.AssemblerScalar;
import edu.uta.futureye.lib.element.FELinearTriangle;
import edu.uta.futureye.util.container.ElementList;
import edu.uta.futureye.util.container.NodeList;
import edu.uta.futureye.util.container.ObjList;

public class LagrangianMethod {
	
	protected double[]s = null;
	protected static String outputFolder = "Lagrangian";
	
	protected Vector u0 = null;
	protected Vector[] g = null;
	protected Vector u0_x = null;
	protected Vector u0_y = null;
	protected Vector lnu0_x = null;
	protected Vector lnu0_y = null;
	protected Vector laplace_ln_u0 = null;
	protected Mesh mesh = null;
	protected Mesh meshBig = null;
	protected Vector a = null; //a(x) from GCM
	protected Vector a_bak = null; //backup of a(x) from GCM
	protected double k;
	protected double theta;
	
	//for test
	protected Vector vReal = null;
	
	/**
	 * 
	 * @param meshBig: big mesh for computing background 'u0'
	 * @param mesh: small mesh
	 * @param s: light sources
	 * @param a: a_glob(x), GCM result
	 * @param k: background of a(x)=k^2
	 * @param g: measurement data (uiSmall or ui)
	 * @param theta: regular param [0,1] or >1?
	 */
	public LagrangianMethod(Mesh meshBig, Mesh mesh, 
			double[]s, Vector a, double k,
			Vector[] g, double theta) {
		this.meshBig = meshBig;
		this.mesh = mesh;
		this.s = s;
		this.a = a;
		this.a_bak = a;
		this.k = k;
		this.g = g;
		this.theta = theta;
		
	}
	
	public void run() {
		ObjList<Vector> vList = new ObjList<Vector>();
		ObjList<Vector> lambdaList = new ObjList<Vector>();
		for(int iter = 0; iter<13; iter++) {
			vList.clear();
			lambdaList.clear();
			for(int i=0;i<s.length;i++) {
				compute_u0(i);
				compute_u0_xy(i);
				compute_Laplace_ln_u0(i);
				Vector v = solve_v(i);
				vList.add(v);
				//v如果是精确解，lambda=0
				Vector lambda = solve_lambda(i, v);
				lambdaList.add(lambda);
			}
			Vector newa = solve_a(iter,vList,lambdaList);
			plotVector(mesh, newa, String.format("Lagrangian_rlt_a%02d.dat",iter));
			a = newa;
			double err = SparseVector.axpy(-1.0, a, a_bak).norm2();
			System.out.println("Iter"+iter+" Error:  Norm2(a(x)-a_glob(x)) =========> "+err);
			br.println("Iter"+iter+" Error:  Norm2(a(x)-a_glob(x)) = "+err);
		}

	}
	
	public Vector solve_a(int iter,ObjList<Vector> vList,ObjList<Vector> lambdaList) {
		double h = s[0] - s[1];
		Vector int_a_b = new SparseVector(vList.at(1).getDim());
		ObjList<Vector> vml = new ObjList<Vector>();
		for(int i=0;i<s.length;i++) {
			//vml_i = v_i*lambda_i
			vml.add(SparseVector.axMuly(1.0, vList.at(i+1), lambdaList.at(i+1)));
			plotVector(mesh, vml.at(i+1), "Lagrangian_vml_"+i+".dat");
		}
		
		for(int i=1;i<s.length;i++) {
			//Trapezoid formula: (a+b)*h/2
			int_a_b = SparseVector.axpy(1.0, int_a_b, 
					SparseVector.ax(h/(2.0*theta),
					SparseVector.axpy(1.0, vml.at(i), vml.at(i+1)))
					);
			plotVector(mesh, int_a_b, "int_a_b_"+i+".dat");
		}
		
		//return int_a_b;
		SparseVector rlt = new SparseVector(this.a.getDim());
		NodeList nodes = mesh.getNodeList();
		for(int i=1;i<=nodes.size();i++) {
			Node node = nodes.at(i);
			double x = node.coord(1);
			double y = node.coord(2);
			//(2.0, 2.6) R=0.3
			if(Math.sqrt((x-2.0)*(x-2.0) + (y-2.6)*(y-2.6)) < 0.15)
				rlt.add(node.globalIndex, this.a.get(node.globalIndex)+int_a_b.get(node.globalIndex));
			else
				rlt.set(node.globalIndex, this.a.get(node.globalIndex));
		}
		return rlt;
		//return SparseVector.axpy(1.0, int_a_b, this.a);
	}
	
	public void compute_u0(int s_i) {
		double h = s[0] - s[1];
		Model model = new Model();
		model.setDelta(1.5+s_i*h, 3.5);
		//Solve background forward problem u0
		model.setMu_a(0.0, 0.0, 0.0, 
				0.1, //mu_a=0.1 mu_s(=model.k)=0.02 => a(x)=5
				1);
		//在mesh上求解，top边界条件不好确定，因此在比较大的meshBig上求解，
		//然后截取出来到较小的mesh上
		Vector u0Big = model.solveForwardNeumann(meshBig);
		plotVector(meshBig, u0Big, "u0_big"+s_i+".dat");
		//截取meshBig的部分解到mesh上
		u0 = model.extractData(meshBig, mesh, u0Big);
		plotVector(mesh, u0, "u0_"+s_i+".dat");
		
		//直接在小网格上求解，自然边界，两者在边界处应该一致，
		//但是经测试不一致，原因是delta函数如果在区域外，结果相差很大
		//Vector u0Samll = model.solveForwardNeumann(mesh);
		Vector u0Samll = model.solveForwardDirichlet(mesh,new Vector2Function(u0));

		plotVector(mesh, u0Samll, "u0_small"+s_i+".dat");
		
		u0=u0Samll;
		
		vReal = SparseVector.axDivy(1.0, g[s_i], u0);
		plotVector(mesh, 
				vReal, "vReal"+s_i+".dat");
	}
	
	public void compute_u0_xy(int s_i) {
		u0_x = Tools.computeDerivative(mesh, u0, "x");
		u0_y = Tools.computeDerivative(mesh, u0, "y");
		plotVector(mesh, u0_x, "u0_x"+s_i+".dat");
		plotVector(mesh, u0_y, "u0_y"+s_i+".dat");
	}
	
	public void compute_Laplace_ln_u0(int s_i) {
		Vector lnu0 = new SparseVector(u0.getDim());
		for(int i=1;i<=u0.getDim();i++) {
			lnu0.set(i, Math.log(u0.get(i)));
		}
		lnu0_x = Tools.computeDerivative(mesh, lnu0, "x");
		lnu0_y = Tools.computeDerivative(mesh, lnu0, "y");
		plotVector(mesh, lnu0_x, "lnu0_x"+s_i+".dat");
		plotVector(mesh, lnu0_y, "lnu0_y"+s_i+".dat");
		
		Vector lnu0_xx = Tools.computeDerivative(mesh, lnu0_x, "x");
		Vector lnu0_yy = Tools.computeDerivative(mesh, lnu0_y, "y");
		laplace_ln_u0 = SparseVector.axpy(1.0, lnu0_xx, lnu0_yy);
		
		plotVector(mesh, laplace_ln_u0, "lnu0_laplace"+s_i+".dat");
	}

	/**
	 * v=w/w0
	 * 
	 * Solve for
	 * \Delta{v} + 2\Nabla{lnu_0}\dot\Nabla{v} - (a-k^2)*v = 0, on \Omega
	 * =>
	 * -\Delta{v} - 2\Nabla{lnu_0}\dot\Nabla{v} + (a-k^2)*v = 0, on \Omega
	 * Robin:
	 * \partial_{n}{v} + (1+\partial_{n}{lnu_0})*v = 1+\partial_{n}{lnu_0}, on \partial\Omega
	 * 
	 * @return
	 */
	public Vector solve_v(int s_i) {
		Vector2Function fu0 = new Vector2Function(u0);
		Vector2Function fu0_x = new Vector2Function(u0_x);
		Vector2Function fu0_y = new Vector2Function(u0_y);
		
		//-( a(x)-k^2 ) = k^2-a(x)
		Vector v_c = SparseVector.axpy(-1.0, a, new SparseVector(a.getDim(),k*k));
		plotVector(mesh, v_c, "param_c"+s_i+".dat");
		Vector2Function param_c = new Vector2Function(v_c);
		
		//\Nabla{lnu0}
		//Function b1 = FOBasic.Divi(fu0_x,fu0);
		//Function b2 = FOBasic.Divi(fu0_y,fu0);
		Function b1 = new Vector2Function(lnu0_x);
		Function b2 = new Vector2Function(lnu0_y);
		plotFunction(mesh, b1, "b1_"+s_i+".dat");
		plotFunction(mesh, b2, "b2_"+s_i+".dat");
		
		WeakFormGCM weakForm = new WeakFormGCM();
		weakForm.setF(FC.c(0.0));
		weakForm.setParam(
				FC.c(-1.0), //注意，有负号!!!
				param_c, 
				FMath.Mult(FC.c(2.0),b1),
				FMath.Mult(FC.c(2.0),b2));
		//q = d = 1 + \partial_{n}{lnu_0}
		Function param_qd = FMath.Plus(
				FC.c(1.0), new DuDn(b1, b2, null)
				);
		//Robin:  d*u + k*u_n = q (自然边界：d==k, q=0)
		weakForm.setRobin(FMath.Mult(FC.c(-1.0),param_qd),
				FMath.Mult(FC.c(-1.0),param_qd));

		mesh.clearBorderNodeMark();
		HashMap<NodeType, Function> mapNTF = new HashMap<NodeType, Function>();
		mapNTF.put(NodeType.Robin, null);
		mesh.markBorderNode(mapNTF);

		AssemblerScalar assembler = new AssemblerScalar(mesh, weakForm);
		System.out.println("Begin Assemble...solve_v");
		assembler.assemble();
		Matrix stiff = assembler.getStiffnessMatrix();
		Vector load = assembler.getLoadVector();
		System.out.println("Assemble done!");

		Solver solver = new Solver();
		Vector v = solver.solve(stiff, load);
		plotVector(mesh, v, "Lagrangian_v"+s_i+".dat");
		return v;
	}
	

	/**
	 * Solve for la=lambda
	 *   \Delta{la} - 2\Nabla{lnu_0}\dot\Nabla{la} - ( 2*Laplace(lnu0)+(a-k^2) )*la = 0, on \Omega
	 * =>
	 *   -\Delta{la} + 2\Nabla{lnu_0}\dot\Nabla{la} + ( 2*Laplace(lnu0)+(a-k^2) )*la = 0, on \Omega
	 * Robin:
	 *   \partial_{n}{la} + (1-\partial_{n}{lnu_0})*la = v - g_tidle, on \partial\Omega
	 * 
	 * where g_tidle = g/u0 = ui/u0|_\Gamma
	 * 
	 * 
	 * March. 7th 2011
	 * //测试将dleta函数变得平缓
		delta = new FDelta(this.lightSource,0.1,2e2);
		
		
	 * @param s_i
	 * @param v
	 * @return
	 */
	public Vector solve_lambda(int s_i, Vector v) {
		Vector2Function fu0 = new Vector2Function(u0);
		Vector2Function fu0_x = new Vector2Function(u0_x);
		Vector2Function fu0_y = new Vector2Function(u0_y);
		
		//-2*Laplace(lnu0) - ( a(x)-k^2 ) = -2*Laplace(lnu0) + ( k^2-a(x) )
		Vector2Function param_c = new Vector2Function(
				SparseVector.axpy(-2.0, laplace_ln_u0,
				SparseVector.ax(1.0,SparseVector.axpy(-1, a, 
						new SparseVector(a.getDim(),k*k)))));
		
		//\Nabla{lnu0}
		Function b1 = FMath.Divi(fu0_x,fu0);
		Function b2 = FMath.Divi(fu0_y,fu0);
		WeakFormGCM weakForm = new WeakFormGCM();

		//(v - g)_\partial{\Omega} 
		Vector v_g = SparseVector.axpy(-1.0, SparseVector.axDivy(1.0,g[s_i],u0), v);
		NodeList nodes = mesh.getNodeList();
		//Test: v_g在整个区域上都已知
//		for(int i=1;i<=nodes.size();i++) {
//			Node node = nodes.at(i);
//			if(Math.abs(node.coord(2)-1.0) < 0.2 && Math.abs(node.coord(1)-2) < 0.2)
//				v_g.set(node.globalIndex, 1.0);
//			else
//				v_g.set(node.globalIndex,0.0);
//		}
		
		Vector v_g2 = SparseVector.axpy(-1.0, SparseVector.axDivy(1.0,g[s_i],u0), v);
		//只保留上边界数据，计算范数
		for(int i=1;i<=nodes.size();i++) {
			Node node = nodes.at(i);
			if(Math.abs(node.coord(2)-3.0) > 0.1)
				v_g2.set(node.globalIndex, 0.0);
		}
		System.out.println("v-g on border norm -------------> "+v_g2.norm2());
		br.println("v-g on border norm = "+v_g2.norm2());
		Function fv_g = new Vector2Function(v_g);

		weakForm.setF(FC.c(0.0));
		//Test: v_g在整个区域上都已知
		//weakForm.setF(fv_g.mult(FC.c(-1.0)));//.mult(FC.c(-20.0))///////////////
		
		weakForm.setParam(
				FC.c(-1.0),//注意，有负号!!!
				param_c, 
				FMath.Mult(FC.c(-2.0),b1),
				FMath.Mult(FC.c(-2.0),b2)
			);
		
		
		plotFunction(mesh, fv_g, "Lagrangian_v_g"+s_i+".dat");
		Function param_d = FMath.Minus(
				FC.c(1.0), new DuDn(b1, b2, null)
				);
		
		//q=v-g_tidle, d=\partial_{n}{lnu_0}
		//Robin:  d*u + k*u_n= q (自然边界：d==k, q=0)
		weakForm.setRobin(FMath.Mult(FC.c(-1.0),fv_g),
				FMath.Mult(FC.c(-1.0), param_d));
		//Test: v_g在整个区域上都已知
		//weakForm.setRobin(FC.c(0.0),
		//		FMath.Mult(FC.c(-1.0), param_d));
		
		mesh.clearBorderNodeMark();
		HashMap<NodeType, Function> mapNTF = new HashMap<NodeType, Function>();
		mapNTF.put(NodeType.Robin, null);
		mesh.markBorderNode(mapNTF);

		AssemblerScalar assembler = new AssemblerScalar(mesh, weakForm);
		System.out.println("Begin Assemble...solve_lambda");
		assembler.assemble();
		Matrix stiff = assembler.getStiffnessMatrix();
		Vector load = assembler.getLoadVector();
		System.out.println("Assemble done!");

		Solver solver = new Solver();
		Vector lambda = solver.solve(stiff, load);
		plotVector(mesh, lambda, "Lagrangian_lambda"+s_i+".dat");
		return lambda;
	}
	
	public static void plotVector(Mesh mesh, Vector v, String fileName) {
	    MeshWriter writer = new MeshWriter(mesh);
	    if(!outputFolder.isEmpty()) {
		    File file = new File("./"+outputFolder);
			if(!file.exists()) {
				file.mkdirs();
			}
	    }
	    writer.writeTechplot("./"+outputFolder+"/"+fileName, v);
	}
		
	public static void plotFunction(Mesh mesh, Function fun, String fileName) {
	    NodeList list = mesh.getNodeList();
	    int nNode = list.size();
		Variable var = new Variable();
		Vector v = new SparseVector(nNode);
	    for(int i=1;i<=nNode;i++) {
	    	Node node = list.at(i);
	    	var.setIndex(node.globalIndex);
	    	var.set("x", node.coord(1));
	    	var.set("y", node.coord(2));
	    	v.set(i, fun.value(var));
	    }
	    plotVector(mesh,v,fileName);
	}
	
	public static void test(String gridFileBig, String gridFileSmall) {
		MeshReader readerForward = new MeshReader(gridFileBig);
		Mesh meshBig = readerForward.read2DMesh();
		MeshReader readerGCM = new MeshReader(gridFileSmall);
		Mesh meshSmall = readerGCM.read2DMesh();
		
		//Use element library to assign degree of freedom (DOF) to element
		ElementList eList = meshBig.getElementList();
		FELinearTriangle linearTriangle = new FELinearTriangle();
		for(int i=1;i<=eList.size();i++)
			linearTriangle.assign(eList.at(i));
		meshBig.computeNodeBelongsToElements();
		meshBig.computeNeighborNodes();
		
		eList = meshSmall.getElementList();
		for(int i=1;i<=eList.size();i++)
			linearTriangle.assign(eList.at(i));
		meshSmall.computeNodeBelongsToElements();
		meshSmall.computeNeighborNodes();


		Model model = new Model();
		//Simulated 'a_glob(x)'
		model.setMu_a(2.0, 2.6, 0.3, 
				0.2, //init 'a_glob(x)'
				1);
		Vector a_glob = model.discreteFunction(meshSmall,
				FMath.Divi(model.mu_a,model.k)
				);
		plotVector(meshSmall, a_glob, "a_glob.dat");

		//Number of moving light sources
		int N=10; 
		double[]s = new double[N];
		double h = 0.2;
		s[0] = 3.0;
		for(int i=1; i<N; i++)
			s[i] = s[0] - i*h; //由远到近，s递减
		
		//Simulated border measurement data 'g'
		Vector[] uiBig = new Vector[N]; //在大网格上求解，截取需要的边界
		Vector[] uiSmall = new Vector[N]; //直接在小网格上求解，自然边界，两者应该一致
		Vector[] ui = new Vector[N];
		
		//实际上，应该来自从tail重构出来的a(x)
		model.setMu_a(2.0, 2.6, 0.3, 
				0.4, //max 'mu_a'
				1);
		plotFunction(meshBig,model.mu_a.D(model.k), "a_real.dat");
		for(int i=0;i<N;i++) {
			model.setDelta(1.5+i*h, 3.5);
			plotFunction(meshBig,model.delta, "delta_"+i+".dat");
			
			//在meshGCM上求解，top边界条件不好确定，因此在比较大的meshForward上求解，
			//然后截取出来到较小的meshGCM上
			uiBig[i] = model.solveForwardNeumann(meshBig); 
			plotVector(meshBig, uiBig[i], "ui_big"+i+".dat");
			//截取meshForward的部分解到meshGCM上
			ui[i] = model.extractData(meshBig, meshSmall, uiBig[i]);
			plotVector(meshSmall, ui[i], "ui(=g)"+i+".dat");
			
			//直接在小网格上求解，自然边界，两者在边界处应该一致
			//但是经测试不一致，原因是delta函数如果在区域外，结果相差很大
			//uiSmall[i] = model.solveForwardNeumann(meshSmall);
			uiSmall[i] = model.solveForwardDirichlet(meshSmall,new Vector2Function(ui[i]));
			
			plotVector(meshSmall, uiSmall[i], "ui_small"+i+".dat");
		}
		
		LagrangianMethod lm = new LagrangianMethod(
				meshBig, //big mesh for computing background 'u0'
				meshSmall, //mesh
				s, //light sources
				a_glob, //a_glob(x)
				Math.sqrt(0.1/model.k.value(null)), //background of a(x)=k^2
				uiSmall, //g(uiSmall or ui), g_tidle=g/u0
				0.1 //theta
				);
		lm.beginLog();
		lm.run();
		lm.endLog();
	}
	
	FileOutputStream out = null;
	PrintWriter br = null;
	public void beginLog() {
		try {
			File file = new File(".\\"+outputFolder+"\\output.txt");
			out = new FileOutputStream(file);
			OutputStreamWriter writer = new OutputStreamWriter(out, "UTF-8");
			br = new PrintWriter(writer);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	public void endLog() {
		try {
			if(br != null)
				br.close();
			if(out != null)
				out.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public static void main(String[] args) {
		test("prostate_test3_ex.grd","prostate_test3.grd");
	}
}
