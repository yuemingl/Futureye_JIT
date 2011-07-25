package edu.uta.futureye.application;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import edu.uta.futureye.algebra.SchurComplementLagrangianSolver;
import edu.uta.futureye.algebra.Solver;
import edu.uta.futureye.algebra.SolverJBLAS;
import edu.uta.futureye.algebra.SparseBlockMatrix;
import edu.uta.futureye.algebra.SparseBlockVector;
import edu.uta.futureye.algebra.SparseMatrix;
import edu.uta.futureye.algebra.SparseVector;
import edu.uta.futureye.algebra.intf.BlockMatrix;
import edu.uta.futureye.algebra.intf.BlockVector;
import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.core.DOF;
import edu.uta.futureye.core.DOFOrder;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.core.NodeRefined;
import edu.uta.futureye.core.NodeType;
import edu.uta.futureye.core.Refiner;
import edu.uta.futureye.core.geometry.GeoEntity;
import edu.uta.futureye.function.AbstractFunction;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.DuDn;
import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.basic.Vector2Function;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.function.intf.ScalarShapeFunction;
import edu.uta.futureye.function.operator.FMath;
import edu.uta.futureye.io.MeshReader;
import edu.uta.futureye.io.MeshWriter;
import edu.uta.futureye.lib.assembler.AssemblerScalar;
import edu.uta.futureye.lib.element.FEBilinearRectangle;
import edu.uta.futureye.lib.element.FELinearTriangle;
import edu.uta.futureye.lib.shapefun.SFBilinearLocal2D;
import edu.uta.futureye.lib.shapefun.SFLinearLocal2D;
import edu.uta.futureye.lib.weakform.WeakFormL22D;
import edu.uta.futureye.lib.weakform.WeakFormLaplace2D;
import edu.uta.futureye.util.FutureyeException;
import edu.uta.futureye.util.PairDoubleInteger;
import edu.uta.futureye.util.Utils;
import edu.uta.futureye.util.container.DOFList;
import edu.uta.futureye.util.container.ElementList;
import edu.uta.futureye.util.container.NodeList;

/**
 * Implementation of the paper: 
 *   'A Framework For The Adaptive Finite Element Solution Of Large-Scale Inverse Problems'
 * 
 * Lagrange Multiplier Method
 * 
 * v=u/u0
 * 
 * @author liuyueming
 *
 */
public class VariationGaussNewtonDOT {
	public Mesh mesh;
	public Mesh meshBig;
	
	
	protected static String outputFolder = "Lagrangian_Klibanov";
	public boolean debug = false;

	//Back ground model
	ModelDOT modelBk = new ModelDOT();
	
	//Real inclusion model
	ModelDOT modelReal = new ModelDOT();

	//Real inclusion model
	ModelDOT modelGuess = new ModelDOT();


    Function diri = null;
    
    //对应每个光源s_i的参数，包括测量数据
    public static class ParamOfLightSource {
    	int s_i;
    	//背景数据
    	Vector u0;
    	Vector lnu0;
    	Vector lnu0_x;
    	Vector lnu0_y;
    	Vector laplace_lnu0;
    	//有包含物的数据
    	Vector u; //光强测量数据（只有边界）
    	Vector g_tidle; //= v|_\Gamma = (u/u0)|_\Gamma = g/u0|_\Gamma
    }
    
    //3*mu_a*mu_s' = 
    // background: 3*0.1*50
    // inclusion:  3*1.0*50
    Vector aGlob = null;//GCM方法得到的a_glob(x)
    double model_k = 0.02;//=1.0/mu_s'
    double aBkSquare = 0.1/model_k;//back ground of a(x)
    
    //正则化参数
    double beta = 10;//for \Delta^-2(a-\overline{a})
    //double beta = 0.1;
    

    int iterNum = 0;
    
    //光源x坐标位置数组
	double[] LS;

    /**
     * 构造函数，设定包含物位置
     */
    public VariationGaussNewtonDOT() {
		//Number of moving light sources
		int N=3; 
		LS = new double[N];
		double h = 1.0;
		LS[0] = 2.0;
		for(int i=1; i<N; i++)
			LS[i] = LS[0] + i*h;
		
		//背景mu_a
		modelBk.setMu_a(0.0, 0.0, 0.0, 
				0.1, //mu_a=0.1 mu_s(=model.k)=0.02 => a(x)=5
				1);
//		//有包含物mu_a，真实模型
//		modelReal.setMu_a(2.57, 2.69, 0.3,
//				2.0, //peak value of mu_a
//				1); //Number of inclusions
//		//有包含物mu_a，猜测模型
//		modelGuess.setMu_a(2.57, 2.69, 0.3,
//				1.0, //peak value of mu_a
//				1); //Number of inclusions
		
//		//有包含物mu_a，真实模型
//		modelReal.setMu_a(3.48, 2.70, 0.3,
//				2.0, //peak value of mu_a
//				1); //Number of inclusions
//		//有包含物mu_a，猜测模型
//		modelGuess.setMu_a(3.48, 2.70, 0.3,
//				1.0, //peak value of mu_a
//				1); //Number of inclusions
		
//		//有包含物mu_a，真实模型
//		modelReal.setMu_a(2.39, 2.70, 0.3,
//				2.0, //peak value of mu_a
//				1); //Number of inclusions
//		//有包含物mu_a，猜测模型
//		modelGuess.setMu_a(2.39, 2.70, 0.3,
//				1.0, //peak value of mu_a
//				1); //Number of inclusions
		
		//有包含物mu_a，真实模型
		modelReal.setMu_a(3.0, 2.50, 0.5,
				2.0, //peak value of mu_a
				1); //Number of inclusions
		//有包含物mu_a，猜测模型
		modelGuess.setMu_a(3.0, 2.50, 0.5,
				1.0, //peak value of mu_a
				1); //Number of inclusions

    }
    
    /**
     * 
     * @param s_i：Light source No. (0,1,2...)
     */
    public void reinitModelLight(int s_i) {
		modelBk.setDelta(LS[s_i], 3.5);
		modelBk.lightNum = s_i;
		modelReal.setDelta(LS[s_i], 3.5);
		modelReal.lightNum = s_i;
		modelGuess.setDelta(LS[s_i], 3.5);
		modelGuess.lightNum = s_i;
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
	

	public void readMeshTriangle(){
		String gridFileBig = "prostate_test3_ex.grd";
		String gridFileSmall = "prostate_test3.grd";

        MeshReader readerBig = new MeshReader(gridFileBig);
        MeshReader readerSmall = new MeshReader(gridFileSmall);
        meshBig = readerBig.read2DMesh();
        mesh = readerSmall.read2DMesh();
        meshBig.computeNodeBelongsToElements();
        mesh.computeNodeBelongsToElements();
        mesh.computeNeighborNodes();

        //2.Mark border types
        HashMap<NodeType, Function> mapNTF =
                new HashMap<NodeType, Function>();
        mapNTF.put(NodeType.Dirichlet, null);
        mesh.markBorderNode(mapNTF);
        
        //3.Use element library to assign degrees of
        //  freedom (DOF) to element
        ElementList eList = mesh.getElementList();
//        for(int i=1;i<=eList.size();i++) {
//        	System.out.println(eList.at(i));
//        }
        FELinearTriangle feLT = new FELinearTriangle();
        for(int i=1;i<=eList.size();i++)
            feLT.assignTo(eList.at(i));
  
        ElementList eListBig = meshBig.getElementList();
		for(int i=1;i<=eListBig.size();i++)
			feLT.assignTo(eListBig.at(i));

	}
	
	public void readMeshRectangle(){
		String gridFileBig = "prostate_test7_ex.grd";
		String gridFileSmall = "prostate_test7.grd";

        MeshReader readerBig = new MeshReader(gridFileBig);
        MeshReader readerSmall = new MeshReader(gridFileSmall);
        meshBig = readerBig.read2DMesh();
        mesh = readerSmall.read2DMesh();
        meshBig.computeNodeBelongsToElements();
        mesh.computeNodeBelongsToElements();
        mesh.computeNeighborNodes();
      
        //2.Mark border types
        HashMap<NodeType, Function> mapNTF =
                new HashMap<NodeType, Function>();
        mapNTF.put(NodeType.Dirichlet, null);
        mesh.markBorderNode(mapNTF);
        
        //3.Use element library to assign degrees of
        //  freedom (DOF) to element
        ElementList eList = mesh.getElementList();
        FEBilinearRectangle feLT = new FEBilinearRectangle();
        for(int i=1;i<=eList.size();i++)
            feLT.assignTo(eList.at(i));
  
        ElementList eListBig = meshBig.getElementList();
		for(int i=1;i<=eListBig.size();i++)
			feLT.assignTo(eListBig.at(i));

	}	
	
	public void refineMesh(ElementList eToRefine) {
        meshBig.computeNodeBelongsToElements();
        meshBig.computeNeighborNodes();
        meshBig.computeGlobalEdge();
        meshBig.computeNeighborElements();
        mesh.computeNodeBelongsToElements();
        mesh.computeNeighborNodes();
		mesh.computeGlobalEdge();
		mesh.computeNeighborElements();
		
		ElementList eToRefineBig = new ElementList();
		for(int i=1;i<=eToRefine.size();i++) {
			Element e = meshBig.getElementByNodes(eToRefine.at(i).nodes);
			eToRefineBig.add(e);
		}
		
		System.out.println("Before refine meshBig: Element="+meshBig.getElementList().size()+", Node="+mesh.getNodeList().size());
		Refiner.refineOnce(meshBig, eToRefineBig);
		System.out.println("After refine: Element="+meshBig.getElementList().size()+", Node="+mesh.getNodeList().size());

		System.out.println("Before refine mesh: Element="+mesh.getElementList().size()+", Node="+mesh.getNodeList().size());
		Refiner.refineOnce(mesh, eToRefine);
		System.out.println("After refine: Element="+mesh.getElementList().size()+", Node="+mesh.getNodeList().size());
		
		assignLinearShapFunction(meshBig);
		assignLinearShapFunction(mesh);
	}
	
	/**
	 * Solve background forward problem u0
	 * 
	 * @param s_i
	 * @return
	 */
	public Vector solveU0(int s_i) {
		//在mesh上求解，top边界条件不好确定，因此在比较大的meshBig上求解，
		//然后截取出来到较小的mesh上
		//Vector u0Neumann = modelBk.solveNeumann(mesh);
		//plotVector(mesh, u0Neumann, "u0_neumann"+s_i+".dat");
		
		Vector u0Big = modelBk.solveNeumann(meshBig);
		plotVector(meshBig, u0Big, "u0_big"+s_i+".dat");
		
		//截取meshBig的部分解到mesh上
		Vector u0Samll = Tools.extractData(meshBig, mesh, u0Big);
		plotVector(mesh, u0Samll, "u0_"+s_i+".dat");
		
		//直接在小网格上求解，自然边界，两者在边界处应该一致，
		//但是经测试不一致，原因是delta函数如果在区域外，结果相差很大
		//u0Samll = model.solveForwardNeumann(mesh);
		
		return u0Samll;
	}
	
	
	
	/**
	 * Get stiff matrix and load vector for equation of lambda
	 * 
	 *   \Delta{lmd} - 2\nabla{lnu_0}\dot\nabla{lmd} - ( 2*Laplace(lnu0)+(a-k^2) )*lmd = 0, on \Omega
	 * =>
	 *   -\Delta{lmd} + 2\nabla{lnu_0}\dot\nabla{lmd} + ( 2*Laplace(lnu0)+(a-k^2) )*lmd = 0, on \Omega
	 * 
	 * Boundary Condition: Robin
	 *   \partial_{n}{lmd} + (1-\partial_{n}{lnu_0})*lmd = v - g_tidle, on \partial\Omega
	 * 
	 * where g_tidle = v|_\Gamma = (ui/u0)|_\Gamma = g/u0|_\Gamma
	 * 
	 * @param s_i
	 * @param v
	 * @return
	 */
	public Equation getEqnLambda(int s_i,  
			Vector lnu0_x, Vector lnu0_y, Vector laplace_lnu0,
			Vector a, Vector v, Vector g_tidle) {
		//-2*Laplace(lnu0) - ( a(x)-k^2 ) = -2*Laplace(lnu0) + ( k^2-a(x) )
		Vector2Function param_c = new Vector2Function(
				FMath.axpy(-2.0, laplace_lnu0,
				FMath.ax(1.0,FMath.axpy(-1, a, 
						new SparseVector(a.getDim(),aBkSquare)))));
		
		//\nabla{lnu0}
		Function b1 = new Vector2Function(lnu0_x);
		Function b2 = new Vector2Function(lnu0_y);
		WeakFormGCM weakForm = new WeakFormGCM();

		//(v - g)_\partial{\Omega} 
		Vector v_g = FMath.axpy(-1.0, g_tidle, v);
		//System.out.println("v-g norm -------------> "+v_g.norm2());
		Function fv_g = new Vector2Function(v_g);
		plotFunction(mesh, fv_g, String.format("M%02d_Lagrangian_v_g%02d.dat",s_i,this.iterNum));

		//weakForm.setF(FC.c(0.0));
		//Test: v_g在整个区域上都已知
		weakForm.setF(fv_g.M(FC.c(-1.0)));

		weakForm.setParam(
				FC.c(-1.0),//注意，有负号!!!
				param_c, 
				FC.c(-2.0).M(b1),
				FC.c(-2.0).M(b2)
			);
		
		Function param_d = FC.c(1.0).S(new DuDn(b1, b2, null));
		
		//q=v-g_tidle, d=\partial_{n}{lnu_0}
		//Robin:  d*u + k*u_n= q (自然边界：d==k, q=0)
		//weakForm.setRobin(FC.c(-1.0).M(fv_g),
		//		FC.c(-1.0).M(param_d));
		//Test: v_g在整个区域上都已知
		weakForm.setRobin(FC.c(0.0),
				FC.c(-1.0).M(param_d));
		
		mesh.clearBorderNodeMark();
		HashMap<NodeType, Function> mapNTF = new HashMap<NodeType, Function>();
		//mapNTF.put(NodeType.Robin, null);
		//Test: v_g在整个区域上都已知
		mapNTF.put(NodeType.Dirichlet, null);
		mesh.markBorderNode(mapNTF);

		AssemblerScalar assembler = new AssemblerScalar(mesh, weakForm);
		System.out.println("Begin Assemble...lambda");
		assembler.assemble();
		Matrix stiff = assembler.getStiffnessMatrix();
		Vector load = assembler.getLoadVector();
		assembler.imposeDirichletCondition(FC.c0);
		System.out.println("Assemble done!");

		Equation eqn = new Equation();
		eqn.A = stiff;
		//注意：这是按照Klibanov的公式来的，与Wolfgang的在lambda上差一个负号
		eqn.f = load.scale(-1.0);
		
		return eqn;
	}	
	
	/**
	 * Get stiff matrix and load vector for equation of v
	 * 
	 *   \Delta{v} + 2\nabla{lnu_0}\dot\nabla{v} - (a-k^2)*v = 0, on \Omega
	 * =>
	 *  -\Delta{v} - 2\nabla{lnu_0}\dot\nabla{v} + (a-k^2)*v = 0, on \Omega
	 *  
	 * Boundary Condition: Robin
	 *   \partial_{n}{v} + (1+\partial_{n}{lnu_0})*v = 1+\partial_{n}{lnu_0}, on \partial\Omega
	 * 
	 * where 
	 *   v=u/u0
	 * 
	 * @param lnu0_x
	 * @param lnu0_y
	 * @param a
	 * @param g_tidle: null->Robin
	 * @return
	 */
	public Equation getEqnV(Vector lnu0_x, Vector lnu0_y, 
			Vector a, Vector g_tidle) {
		//-( a(x)-k^2 ) = k^2-a(x)
		Vector v_c = FMath.axpy(-1.0, a, new SparseVector(a.getDim(),aBkSquare));
		//plotVector(mesh, v_c, "param_c"+s_i+".dat");
		Vector2Function param_c = new Vector2Function(v_c);
		
		Function b1 = new Vector2Function(lnu0_x);
		Function b2 = new Vector2Function(lnu0_y);
		
		WeakFormGCM weakForm = new WeakFormGCM();
		weakForm.setF(FC.c(0.0));
		weakForm.setParam(
				FC.c(-1.0), //注意，有负号!!!
				param_c, 
				FC.c(2.0).M(b1),
				FC.c(2.0).M(b2));
		
		//q = d = 1 + \partial_{n}{lnu_0}
		Function param_qd = FC.c(1.0).A(new DuDn(b1, b2, null));
		weakForm.setRobin(FC.c(-1.0).M(param_qd),
				FC.c(-1.0).M(param_qd));

		mesh.clearBorderNodeMark();
		HashMap<NodeType, Function> mapNTF = new HashMap<NodeType, Function>();
		if(g_tidle == null)
			mapNTF.put(NodeType.Robin, null);
		else
			mapNTF.put(NodeType.Dirichlet, null);
		mesh.markBorderNode(mapNTF);

		AssemblerScalar assembler = new AssemblerScalar(mesh, weakForm);
		System.out.println("Begin Assemble...v");
		assembler.assemble();
		Matrix stiff = assembler.getStiffnessMatrix();
		Vector load = assembler.getLoadVector();
		if(g_tidle != null)
			assembler.imposeDirichletCondition(new Vector2Function(g_tidle));
		System.out.println("Assemble done!");

		Equation eqn = new Equation();
		eqn.A = stiff;
		eqn.f = load;
		
		return eqn;
	}	
	
	/**
	 * Get stiff matrix and load vector for equation of a(x)
	 * 
	 * \beta(a-a_glob,\chi) = -(\chi*v,\lambda)
	 * =>
	 * (a,\chi) = ( -\frac{1.0}{\beta}*v*\lambda + a_glob,\chi )
	 * 
	 * Boundary Condition: Dirichlet
	 * a(x)|_\Gamma = diri (=0.1 Background)
	 * 
	 * @param v
	 * @param lambda
	 * @param a_glob
	 * @return
	 */
	public Equation getEqnA(Vector v, Vector lambda, Vector a_glob, Function diri) {
        //4.Weak form
        WeakFormL22D weakForm = new WeakFormL22D();
        weakForm.setParam(FC.c0, FC.c1);
        //Right hand side(RHS): f(x) = -\frac{1.0}{\beta}*v*\lambda + a_glob
        Function fv = new Vector2Function(v);
        Function flmd = new Vector2Function(lambda);
        Function fa_glob = new Vector2Function(a_glob);
        Function fv_lmd = fv.M(flmd);
        plotFunction(mesh,fv_lmd,String.format("La_RHS_v_lmd%02d.dat",this.iterNum));
//        Function f2 = fv_lmd.M(-1.0/beta).A(fa_glob);
        Function f2 = fv_lmd.M(1.0/beta).A(fa_glob);
        plotFunction(mesh,f2,String.format("La_RHS_all%02d.dat",this.iterNum));
        weakForm.setF(f2);

        //需要重新标记边界条件
		mesh.clearBorderNodeMark();
		HashMap<NodeType, Function> mapNTF = new HashMap<NodeType, Function>();
		mapNTF.put(NodeType.Dirichlet, null);
		mesh.markBorderNode(mapNTF);
		
        //5.Assembly process
        AssemblerScalar assembler =
                new AssemblerScalar(mesh, weakForm);
        System.out.println("Begin Assemble...a(x)");
        assembler.assemble();
        Matrix stiff = assembler.getStiffnessMatrix();
        Vector load = assembler.getLoadVector();
        //Boundary condition
        assembler.imposeDirichletCondition(diri);
        System.out.println("Assemble done!");

		Equation eqn = new Equation();
		eqn.A = stiff;
		eqn.f = load;
		
		return eqn;
	}
	
	/**
	 * \beta(a-a_glob,\chi) = -\sum_{i=1,N}{(\chi*vi,{\lambda}i)}
	 */
	public Equation getEqnA(Vector[] v, Vector[] lambda, Vector a_glob, Function diri) {
        //4.Weak form
        WeakFormL22D weakForm = new WeakFormL22D();
        weakForm.setParam(FC.c0, FC.c1);//????是否该使用WeakFormLaplace2D
        
        //Right hand side(RHS): f(x) = -\frac{1.0}{\beta}*v*\lambda + a_glob
        Function fa_glob = new Vector2Function(a_glob);
        int N=v.length;
        Vector[] vMlmd = new Vector[N];
        for(int i=0; i<N; i++) {
        	//v[i].axMuly(1.0,lambda[i]) 会改变v[i]的值
        	vMlmd[i] = FMath.axMuly(1.0, v[i],lambda[i]);
        	
        	plotVector(mesh,vMlmd[i],
        			String.format("M%02d_La_RHS_v_lmd%02d.dat",i,this.iterNum));
        }
        Vector sum = FMath.sum(vMlmd);
        plotVector(mesh,sum,String.format("La_RHS_v_lmd_sum%02d.dat",this.iterNum));
        
        ///////////////////////////////////////////////
        //\Delta(\Delta(sum))
        sum = Tools.computeLaplace2D(mesh, sum);
        sum = Tools.computeLaplace2D(mesh, sum);
        ///////////////////////////////////////////////
        
        plotVector(mesh,sum,String.format("La_RHS_v_lmd_sum_LL%02d.dat",this.iterNum));
        //Function f2 = new Vector2Function(sum).M(-1.0/beta).A(fa_glob);
        //Function f2 = new Vector2Function(sum).M(1.0/beta).A(fa_glob);
        //？？？1.0/beta=10???
        Function f2 = new Vector2Function(sum).M(1.0/beta).A(fa_glob);
        
        plotFunction(mesh,f2,String.format("La_RHS_all%02d.dat",this.iterNum));
        weakForm.setF(f2);

        //需要重新标记边界条件
		mesh.clearBorderNodeMark();
		HashMap<NodeType, Function> mapNTF = new HashMap<NodeType, Function>();
		mapNTF.put(NodeType.Dirichlet, null);
		mesh.markBorderNode(mapNTF);
		
        //5.Assembly process
        AssemblerScalar assembler =
                new AssemblerScalar(mesh, weakForm);
        System.out.println("Begin Assemble...a(x)");
        assembler.assemble();
        Matrix stiff = assembler.getStiffnessMatrix();
        Vector load = assembler.getLoadVector();
        //Boundary condition
        assembler.imposeDirichletCondition(diri);
        System.out.println("Assemble done!");

		Equation eqn = new Equation();
		eqn.A = stiff;
		eqn.f = load;
		
		return eqn;
	}
	
	/**
	 * State equation L_{\lambda}
	 * 
	 * 获取状态方程的余量
	 * @return
	 */
	public Vector getResLlmd(Vector lnu0_x, Vector lnu0_y, 
			Vector a, Vector v, Vector g_tidle) {
		Equation eq = this.getEqnV(lnu0_x, lnu0_y, a, g_tidle);
        //Residual
        Vector res = new SparseVector(eq.f.getDim());
        eq.A.mult(v, res);
        res.add(-1.0, eq.f);
        return res;
	}
	
	/**
	 * Adjoint equation L_{u}
	 * 
	 * @return
	 */
	public Vector getResLu(int s_i,  
			Vector lnu0_x, Vector lnu0_y, Vector laplace_lnu0,
			Vector a, Vector v, Vector g_tidle,
			Vector lambda) {
		Equation eq = this.getEqnLambda(s_i, lnu0_x, lnu0_y, laplace_lnu0, a, v, g_tidle);
        //Residual
        Vector res = new SparseVector(eq.f.getDim());
        eq.A.mult(lambda, res);
        res.add(-1.0, eq.f);
        return res;
    }
	
	/**
	 * Parameter regularization L_{q}
	 * 
	 * @param v
	 * @param lambda
	 * @param a_glob
	 * @return
	 */
	public Vector getResLq(Vector a_glob, 
			Vector v, Vector lambda, Vector a, 
			Function diri) {
		Equation eq = this.getEqnA(v, lambda, a_glob, diri);
        //Residual
        Vector res = new SparseVector(eq.f.getDim());
        eq.A.mult(a, res);
        res.add(-1.0, eq.f);
        plotVector(mesh,res,String.format("Res_La%02d.dat",this.iterNum));
        res = Utils.gaussSmooth(mesh, res, 2, 0.5);
        res = Utils.gaussSmooth(mesh, res, 2, 0.5);
        plotVector(mesh,res,String.format("Res_LaSmooth%02d.dat",this.iterNum));
        Solver sol = new Solver();
        Vector a_solve = sol.solveCGS(eq.A, eq.f);
        plotVector(mesh,a_solve,String.format("a_solve%02d.dat",this.iterNum));
        return res;
	}
	
	public Vector getResLq(Vector a_glob, 
			Vector[] v, Vector[] lambda, Vector a, 
			Function diri) {
		Equation eq = this.getEqnA(v, lambda, a_glob, diri);
        //Residual
        Vector res = new SparseVector(eq.f.getDim());
        eq.A.mult(a, res);
        res.add(-1.0, eq.f);
        plotVector(mesh,res,String.format("Res_La%02d.dat",this.iterNum));
        
        //光滑余量
        res = Utils.gaussSmooth(mesh, res, 1, 0.5);
        res = Utils.gaussSmooth(mesh, res, 1, 0.5);
        res = Utils.gaussSmooth(mesh, res, 1, 0.5);
        res = Utils.gaussSmooth(mesh, res, 1, 0.5);
        plotVector(mesh,res,String.format("Res_LaSmooth%02d.dat",this.iterNum));
        
        //直接求解a(x)
        Solver sol = new Solver();
        Vector a_solve = sol.solveCGS(eq.A, eq.f);
        plotVector(mesh,a_solve,String.format("a_solve%02d.dat",this.iterNum));
        return res;
	}
	
	/**
	 * Search direction matrix:
	 * 
	 *  ( M  A'  0 )(du)     (Lu)
	 *  ( A  0   C )(dl) = - (Ll)
	 *  ( 0  C' bR )(dq)     (Lq)
	 * 
	 * @return
	 */
	public Matrix getM() {
        //4.Weak form: (du,\phi)
        WeakFormLaplace2D weakForm = new WeakFormLaplace2D();
        weakForm.setParam(FC.c0, FC.c1, null, null);
        //Right hand side(RHS): f(x) = 0
        weakForm.setF(FC.c0);

        //需要重新标记边界条件，否则在“整体合成”过程会出现错误。
        //虽然边界条件实在大矩阵中设置，这一步也是需要的。
		mesh.clearBorderNodeMark();
		HashMap<NodeType, Function> mapNTF = new HashMap<NodeType, Function>();
		mapNTF.put(NodeType.Dirichlet, null);
		mesh.markBorderNode(mapNTF);

        //5.Assembly process
        AssemblerScalar assembler =
                new AssemblerScalar(mesh, weakForm);
        System.out.println("Begin Assemble...M");
        assembler.assemble();
        Matrix stiff = assembler.getStiffnessMatrix();
        //Boundary condition 
        //边界条件需要在大矩阵中设置
        //assembler.imposeDirichletCondition(FC.c0);
        System.out.println("Assemble done!");
        
        return stiff;
    }
	
	/**
	 * (\nabla{dl},\nabla{\phi}) - (2*\nabla{lnu_0}\cdot\nabla{dl},\phi) + ((a(x)-k^2)dl,\phi)
	 * where 
	 *   dl=\delta\lambda
	 */
	public Matrix getA(Vector ak, Vector lnu0_x, Vector lnu0_y) {
        //4.Weak form:
        WeakFormGCM weakForm = new WeakFormGCM();
        Function fak = new Vector2Function(ak);
        
        Function flnu0_x = new Vector2Function(lnu0_x);
        Function flnu0_y = new Vector2Function(lnu0_y);
        //weakForm.setParam(FC.c(k), fak.S(k*k), flnu0_x.M(-2.0), flnu0_y.M(-2.0));
        weakForm.setParam(FC.c1, fak.S(aBkSquare), flnu0_x.M(-2.0), flnu0_y.M(-2.0));
        
        //Right hand side(RHS): f(x) = 0
        weakForm.setF(FC.c0);

        //需要重新标记边界条件，否则在“整体合成”过程会出现错误。
        //虽然边界条件实在大矩阵中设置，这一步也是需要的。
		mesh.clearBorderNodeMark();
		HashMap<NodeType, Function> mapNTF = new HashMap<NodeType, Function>();
		mapNTF.put(NodeType.Dirichlet, null);
		mesh.markBorderNode(mapNTF);
		
        //5.Assembly process
        AssemblerScalar assembler =
                new AssemblerScalar(mesh, weakForm);
        System.out.println("Begin Assemble...A");
        assembler.assemble();
        Matrix stiff = assembler.getStiffnessMatrix();
        //Boundary condition
        //边界条件需要在大矩阵中设置
        //assembler.imposeDirichletCondition(FC.c0);
        System.out.println("Assemble done!");
        
        return stiff;
	}
	
	/**
	 * 
	 * @param ak
	 * @param lnu0_x
	 * @param lnu0_y
	 * @param _resLlmd_da: -resLl-\delta{a}
	 * @param vk
	 * @return
	 */
	public Vector solveDeltaV(
			Vector ak, Vector lnu0_x, Vector lnu0_y,
			Vector _resLlmd_da, Vector vk
			) {
        //4.Weak form:
        WeakFormGCM weakForm = new WeakFormGCM();
        Function fak = new Vector2Function(ak);
        
        Function flnu0_x = new Vector2Function(lnu0_x);
        Function flnu0_y = new Vector2Function(lnu0_y);
        //weakForm.setParam(FC.c(k), fak.S(k*k), flnu0_x.M(-2.0), flnu0_y.M(-2.0));
        weakForm.setParam(FC.c1, fak.S(aBkSquare), 
        		flnu0_x.M(-2.0), flnu0_y.M(-2.0));
        
        //Right hand side(RHS): f(x) = da*vk
        //Vector dqvk = _resLlmd_da.axMuly(1.0, vk);
        Vector dqvk = FMath.axMuly(1.0, _resLlmd_da, vk);
        weakForm.setF(new Vector2Function(dqvk));

        //需要重新标记边界条件，否则在“整体合成”过程会出现错误。
        //虽然边界条件实在大矩阵中设置，这一步也是需要的。
		mesh.clearBorderNodeMark();
		HashMap<NodeType, Function> mapNTF = new HashMap<NodeType, Function>();
		mapNTF.put(NodeType.Dirichlet, null);
		mesh.markBorderNode(mapNTF);
		
        //5.Assembly process
        AssemblerScalar assembler =
                new AssemblerScalar(mesh, weakForm);
        System.out.println("Begin Assemble...solveDeltaV");
        assembler.assemble();
        Matrix stiff = assembler.getStiffnessMatrix();
        Vector load = assembler.getLoadVector();
        //Boundary condition
        assembler.imposeDirichletCondition(FC.c0);
        System.out.println("Assemble done!");
        
        Solver sol = new Solver();
        Vector x = sol.solveCGS(stiff, load);
        return x;
	}
	
	public Matrix getAT(Vector ak, Vector lnu0_x, Vector lnu0_y) {
        //4.Weak form:
        WeakFormGCMDual weakForm = new WeakFormGCMDual();
        Function fak = new Vector2Function(ak);
        
        Function flnu0_x = new Vector2Function(lnu0_x);
        Function flnu0_y = new Vector2Function(lnu0_y);
        //weakForm.setParam(FC.c(k), fak.S(k*k), flnu0_x.M(-2.0), flnu0_y.M(-2.0));
        weakForm.setParam(FC.c1, fak.S(aBkSquare), flnu0_x.M(-2.0), flnu0_y.M(-2.0));
        
        //Right hand side(RHS): f(x) = 0
        weakForm.setF(FC.c0);

        //需要重新标记边界条件，否则在“整体合成”过程会出现错误。
        //虽然边界条件实在大矩阵中设置，这一步也是需要的。
		mesh.clearBorderNodeMark();
		HashMap<NodeType, Function> mapNTF = new HashMap<NodeType, Function>();
		mapNTF.put(NodeType.Dirichlet, null);
		mesh.markBorderNode(mapNTF);
		
        //5.Assembly process
        AssemblerScalar assembler =
                new AssemblerScalar(mesh, weakForm);
        System.out.println("Begin Assemble...AT");
        assembler.assemble();
        Matrix stiff = assembler.getStiffnessMatrix();
        //Boundary condition
        //边界条件需要在大矩阵中设置
        //assembler.imposeDirichletCondition(FC.c0);
        System.out.println("Assemble done!");
        
        return stiff;
	}
	
	public Matrix getC(Vector vk) {
        //4.Weak form: (\phi,da*vk)
        WeakFormLaplace2D weakForm = new WeakFormLaplace2D();
        Function fvk = new Vector2Function(vk);
        weakForm.setParam(FC.c0, fvk, null, null);
        //Right hand side(RHS): f(x) = 0
        weakForm.setF(FC.c0);

        //需要重新标记边界条件，否则在“整体合成”过程会出现错误。
        //虽然边界条件实在大矩阵中设置，这一步也是需要的。
		mesh.clearBorderNodeMark();
		HashMap<NodeType, Function> mapNTF = new HashMap<NodeType, Function>();
		mapNTF.put(NodeType.Dirichlet, null);
		mesh.markBorderNode(mapNTF);
		
        //5.Assembly process
        AssemblerScalar assembler =
                new AssemblerScalar(mesh, weakForm);
        System.out.println("Begin Assemble...C");
        assembler.assemble();
        Matrix stiff = assembler.getStiffnessMatrix();
        //Boundary condition
        //边界条件需要在大矩阵中设置
        //assembler.imposeDirichletCondition(FC.c0);
        System.out.println("Assemble done!");
        
        return stiff;
	}
	
	public Matrix getBR() {
        //4.Weak form: (dq,\chi)
        WeakFormLaplace2D weakForm = new WeakFormLaplace2D();
        weakForm.setParam(FC.c0, FC.c(beta), null, null);
        //Right hand side(RHS): f(x) = 0
        weakForm.setF(FC.c0);

        //需要重新标记边界条件，否则在“整体合成”过程会出现错误。
        //虽然边界条件实在大矩阵中设置，这一步也是需要的。
		mesh.clearBorderNodeMark();
		HashMap<NodeType, Function> mapNTF = new HashMap<NodeType, Function>();
		mapNTF.put(NodeType.Dirichlet, null);
		mesh.markBorderNode(mapNTF);
		
        //5.Assembly process
        AssemblerScalar assembler =
                new AssemblerScalar(mesh, weakForm);
        System.out.println("Begin Assemble...R");
        assembler.assemble();
        Matrix stiff = assembler.getStiffnessMatrix();
        //Boundary condition
        //边界条件需要在大矩阵中设置
        //assembler.imposeDirichletCondition(FC.c0);
        System.out.println("Assemble done!");
        
        return stiff;
	}
	
	/**
	 * v=u/u0
	 * * u|_\Gamma = g
	 * * v|_\Gamma = g_tidle
	 * @param lnu0_x
	 * @param lnu0_y
	 * 
	 * @return
	 */
	public Vector solveRealV(Vector lnu0_x, Vector lnu0_y) {
		NodeList nodes = mesh.getNodeList();
		Vector aReal = new SparseVector(nodes.size());
		for(int i=1;i<=nodes.size();i++) {
			Node node = nodes.at(i);
			double mu_a = modelReal.mu_a.value(Variable.createFrom(modelReal.mu_a, node, i));
			aReal.set(i,mu_a/this.model_k);
		}
		plotVector(mesh,aReal,"aReal.dat");
		Equation eq = this.getEqnV(lnu0_x, lnu0_y, aReal, null);
        Solver solver = new Solver();
        Vector v = solver.solveCGS(eq.A, eq.f);
        return v;
	}
	
	/**
	 * 求解关于v的状态方程
	 * 
	 * @param lnu0_x
	 * @param lnu0_y
	 * @param a 
	 * @return
	 */
	public Vector solveStateEquation(Vector lnu0_x, Vector lnu0_y, Vector a,
			Vector g_tidle) {
		Equation eq = this.getEqnV(lnu0_x, lnu0_y, a, g_tidle);
        Solver solver = new Solver();
        Vector v = solver.solveCGS(eq.A, eq.f);
        return v;
	}
	
	/**
	 * 求解关于lambda的伴随方程
	 * 
	 * @param s_i
	 * @param lnu0_x
	 * @param lnu0_y
	 * @param laplace_lnu0
	 * @param a
	 * @param v
	 * @param g_tidle: v|_\Gamma = (u/u0)|_\Gamma = g/u0|_\Gamma
	 * @return
	 */
	public Vector solveAdjointEquation(int s_i,  
			Vector lnu0_x, Vector lnu0_y, Vector laplace_lnu0,
			Vector a, Vector v, Vector g_tidle) {
		Equation eq = this.getEqnLambda(s_i, lnu0_x, lnu0_y, laplace_lnu0, a, v, g_tidle);
        Solver solver = new Solver();
        Vector lmd_solve = solver.solveCGS(eq.A, eq.f);
        return lmd_solve;
    }
	
	protected void setDirichlet(BlockMatrix BM, BlockVector BV,
			int matIndex, double value) {
		int row = matIndex;
		int col = matIndex;
		BM.set(row, col, 1.0);
		BV.set(row,value);
		for(int r=1;r<=BM.getRowDim();r++) {
			if(r != row) {
				BV.add(r,-BM.get(r, col)*value);
				BM.set(r, col, 0.0);
			}
		}
		for(int c=1;c<=BM.getColDim();c++) {
			if(c != col) {
				BM.set(row, c, 0.0);
			}
		}
	}
	
	public void imposeDirichletCondition(BlockMatrix BM, BlockVector BV,
			Function diri) {
		ElementList eList = mesh.getElementList();
		int nNode = mesh.getNodeList().size();
		for(int i=1;i<=eList.size();i++) {
			Element e = eList.at(i);
			DOFList DOFs = e.getAllDOFList(DOFOrder.NEFV);
			for(int j=1;j<=DOFs.size();j++) {
				DOF dof = DOFs.at(j);
				GeoEntity ge = dof.getOwner();
				if(ge instanceof Node) {
					Node n = (Node)ge;
					if(n.getNodeType() == NodeType.Dirichlet) {
						Variable v = Variable.createFrom(diri, n, 0);
						setDirichlet(BM,BV,dof.getGlobalIndex(),diri.value(v));
						//setDirichlet(BM,BV,nNode+dof.getGlobalIndex(),diri.value(v));
						//setDirichlet(BM,BV,nNode*2+dof.getGlobalIndex(),diri.value(v));
					}
				}
			}
		}
	}
	
	public void imposeDirichletCondition(BlockMatrix BM, BlockVector BV,
			int nDataBlock, Function diri) {
		ElementList eList = mesh.getElementList();
		int nNode = mesh.getNodeList().size();
		for(int i=1;i<=eList.size();i++) {
			Element e = eList.at(i);
			DOFList DOFs = e.getAllDOFList(DOFOrder.NEFV);
			for(int j=1;j<=DOFs.size();j++) {
				DOF dof = DOFs.at(j);
				GeoEntity ge = dof.getOwner();
				if(ge instanceof Node) {
					Node n = (Node)ge;
					if(n.getNodeType() == NodeType.Dirichlet) {
						Variable v = Variable.createFrom(diri, n, 0);
						for(int k=0;k<nDataBlock;k++) {
							setDirichlet(BM,BV,k*nNode+dof.getGlobalIndex(),diri.value(v));
							setDirichlet(BM,BV,(nDataBlock+k)*nNode+dof.getGlobalIndex(),diri.value(v));
						}
						setDirichlet(BM,BV,nDataBlock*2*nNode+dof.getGlobalIndex(),diri.value(v));
					}
				}
			}
		}
	}
	
	/**
	 * 
	 *  (M  A'  0)
	 *  (A  0   C)
	 *  (0  C'  R)
	 * 
	 */
	public BlockMatrix getSearchDirectionMatrix(Vector ak, Vector vk,
			Vector lnu0_x, Vector lnu0_y) {
		BlockMatrix BM = new SparseBlockMatrix(3,3);
		Matrix M = this.getM();
		//A <--> AT ?? fix
		Matrix A = this.getA(ak, lnu0_x, lnu0_y);
		//Matrix AT = A.copy().trans();
		Matrix AT = this.getAT(ak, lnu0_x, lnu0_y);
		Matrix C = this.getC(vk);
		Matrix CT = C.copy().trans();
		Matrix R = this.getBR();
		
		Matrix BM13 = new SparseMatrix(M.getRowDim(),R.getColDim());
		Matrix BM22 = new SparseMatrix(A.getRowDim(),A.getColDim());
		Matrix BM31 = new SparseMatrix(R.getRowDim(),M.getColDim());
		
		BM.setBlock(1, 1, M);
		BM.setBlock(1, 2, AT);
		BM.setBlock(1, 3, BM13);
		
		BM.setBlock(2, 1, A);
		BM.setBlock(2, 2, BM22);
		BM.setBlock(2, 3, C);
		
		BM.setBlock(3, 1, BM31);
		BM.setBlock(3, 2, CT);
		BM.setBlock(3, 3, R);
		
		return BM;
	}
	
	/**
	 * 
	 * (M1        |AT1          | 0  )
	 * (  M2      |   AT2       | 0  )
	 * (    ...   |      ...    | .. )
	 * (       MN |         ATN | 0  )
	 *  ----------------------------
	 * (A1        |0            | C1 )
	 * (  A2      |    0        | C2 )
	 * (    ...   |      ...    | .. )
	 * (       AN |          0  | CN )
	 *  ----------------------------
	 * (0 0 ... 0 |CT1 CT2 . CTN| R  )
	 * 
	 */	
	public BlockMatrix getSearchDirectionMatrix(Vector ak, Vector[] vk,
			List<ParamOfLightSource> paramList) {
		
		int nDataBlock = vk.length;
		BlockMatrix BM = new SparseBlockMatrix(nDataBlock*2+1,nDataBlock*2+1);
		
		Matrix R = this.getBR();
		BM.setBlock(nDataBlock*2+1, nDataBlock*2+1, R.setName("R"));
		
		for(int i=1;i<=nDataBlock;i++) {
			ParamOfLightSource param=paramList.get(i-1);
			Matrix M = this.getM();
			Matrix A = this.getA(ak, param.lnu0_x, param.lnu0_y);
			//Matrix AT = A.copy().trans();
			Matrix AT = this.getAT(ak, param.lnu0_x, param.lnu0_y);
			Matrix C = this.getC(vk[i-1]);
			Matrix CT = C.copy().trans();
			for(int j=1;j<=nDataBlock;j++) {
				if(i==j) {
					BM.setBlock(i, j, M.setName("M"+i));
					BM.setBlock(i, nDataBlock+j, AT.setName("AT"+i));
					BM.setBlock(nDataBlock+i, j, A.setName("A"+i));
					BM.setBlock(nDataBlock+i, nDataBlock*2+1, C.setName("C"+i));
					BM.setBlock(nDataBlock*2+1, nDataBlock+i, CT.setName("CT"+i));
					
					BM.setBlock(nDataBlock+i, nDataBlock+i, 
							new SparseMatrix(A.getRowDim(),AT.getColDim()));
					BM.setBlock(i, nDataBlock*2+1, 
							new SparseMatrix(M.getRowDim(),R.getColDim()));
					BM.setBlock(nDataBlock*2+1, i, 
							new SparseMatrix(R.getRowDim(),M.getColDim()));
					
				} else {
					
					Matrix M0 = new SparseMatrix(M.getRowDim(),M.getColDim());
					Matrix AT0 = new SparseMatrix(AT.getRowDim(),AT.getColDim());
					Matrix A0 = new SparseMatrix(A.getRowDim(),A.getColDim());
					Matrix ATA0 = new SparseMatrix(A.getRowDim(),AT.getColDim());
					BM.setBlock(i, j, M0);
					BM.setBlock(i, nDataBlock+j, AT0);
					BM.setBlock(nDataBlock+i, j, A0);
					BM.setBlock(nDataBlock+i, nDataBlock+j, ATA0);
				}
			}		
		}
		return BM;
	}
	
	/**
	 * Return a new BolckMatrix object that share the same sub matrix objects
	 * 
	 * @param BM
	 * @param col1
	 * @param col2
	 * @return
	 */
	public BlockMatrix changeBlockColumn(BlockMatrix BM, int col1, int col2) {
		int blockRow = BM.getRowBlockDim();
		int blockCol = BM.getColBlockDim();
		
		BlockMatrix newBM = new SparseBlockMatrix(blockRow, blockCol);
		for(int i=1;i<=blockRow;i++) {
			for(int j=1;j<=blockCol;j++) {
				newBM.setBlock(i, j, BM.getBlock(i, j));
			}
		}
		
		for(int i=1;i<=BM.getRowBlockDim();i++) {
			newBM.setBlock(i, col1, BM.getBlock(i, col2));
			newBM.setBlock(i, col2, BM.getBlock(i, col1));
		}
		return newBM;
	}	
	
	public BlockMatrix changeBlockColumn(BlockMatrix BM, int nDataBlock) {
		int blockRow = BM.getRowBlockDim();
		int blockCol = BM.getColBlockDim();
		
		BlockMatrix newBM = new SparseBlockMatrix(blockRow, blockCol);
		for(int i=1;i<=blockRow;i++) {
			for(int j=1;j<=blockCol;j++) {
				newBM.setBlock(i, j, BM.getBlock(i, j));
			}
		}
		
		for(int i=1;i<=BM.getRowBlockDim();i++) {
			for(int col=1;col<=nDataBlock;col++) {
				newBM.setBlock(i, col, BM.getBlock(i, nDataBlock+col));
				newBM.setBlock(i, nDataBlock+col, BM.getBlock(i, col));
			}
		}
		return newBM;
	}		
	
	/**
	 * 
	 * @param paramList
	 */
	public void gaussNewtonIterate(List<ParamOfLightSource> paramList) {
		if(paramList.size()>1) {
			throw new FutureyeException("Unsupported data!");
		}

		ParamOfLightSource param =  paramList.get(0);
		
		//*************************Initial Values********************
		//a0=a_glob
		Vector a0 = this.aGlob;
		//v0
		Vector v0 = this.solveStateEquation(param.lnu0_x,param.lnu0_y,a0,null);
		plotVector(mesh, v0, "v0.dat");
		//lambda0
		Vector lambda0 = this.solveAdjointEquation(param.s_i, param.lnu0_x, param.lnu0_y, 
				param.laplace_lnu0, a0, v0, param.g_tidle);
		plotVector(mesh, lambda0, "lambda0.dat");
		
		//************************************************************
		
		//迭代求解
		BlockVector f = new SparseBlockVector(3);
		for(int iter=0;iter<10;iter++) {
			
			this.iterNum = iter;
			
			BlockMatrix BM = this.getSearchDirectionMatrix(a0, v0,
					param.lnu0_x, param.lnu0_y);
			
			//状态方程余量
			Vector resLlmd = this.getResLlmd(param.lnu0_x,  param.lnu0_y,  
					a0, v0,param.g_tidle).scale(-1.0);
	        plotVector(mesh,resLlmd,
	        		String.format("Res_Llambda%02d.dat",this.iterNum));
			f.setBlock(2, resLlmd);			
			
			//伴随方程余量
			Vector resLu = this.getResLu(param.s_i, 
				param.lnu0_x, param.lnu0_y, param.laplace_lnu0,
				a0, v0, param.g_tidle,lambda0).scale(-1.0);//？1.0
	        plotVector(mesh,resLu,
	        		String.format("Res_Lv%02d.dat",this.iterNum));
			f.setBlock(1, resLu);

			//正则化参数方程余量
			f.setBlock(3, this.getResLq(aGlob, 
					v0, lambda0, a0, 
					FC.c(0.1) //a(x)的边界条件=背景mu_a
					).scale(-1.0));

			//设置边界条件并求解
			//需要交换矩阵BM的第1, 2列，然后设置边界条件，这样使得矩阵A, A'可逆
			BlockMatrix newBM = this.changeBlockColumn(BM, 1, 2);
			this.imposeDirichletCondition(newBM, f, FC.c0);
			SchurComplementLagrangianSolver solver = 
				new SchurComplementLagrangianSolver(BM, f,mesh);
			BlockVector x = solver.solve();
			
			this.imposeDirichletCondition(BM, f, FC.c0);
			//SolverJBLAS sol = new SolverJBLAS();
			//BlockVector x = (BlockVector)sol.solveDGESV(BM, f);
			
			plotVector(mesh, x.getBlock(1), String.format("delta_v%02d.dat",iter));
			plotVector(mesh, x.getBlock(2), String.format("delta_lambda%02d.dat",iter));
			plotVector(mesh, x.getBlock(3), String.format("delta_a%02d.dat",iter));
			
			//待定，与beta选取有关
			double stepLength = 0.1;
			
			v0.add(stepLength, x.getBlock(1));
			lambda0.add(stepLength, x.getBlock(2));
			a0.add(stepLength, x.getBlock(3));
			
			plotVector(mesh, v0, String.format("rlt_v%02d.dat",iter));
			plotVector(mesh, lambda0, String.format("rlt_lambda%02d.dat",iter));
			plotVector(mesh, a0, String.format("rlt_a%02d.dat",iter));

		}

		
	}
	
	
	public Vector gaussNewtonIterateMulti(List<ParamOfLightSource> paramList) {
		int nDataBlock = paramList.size();
		
		//*************************Initial Values********************
		//a0=a_glob
		//copy很重要
		Vector a0 = this.aGlob.copy();
//		Function fa0 = new AbstractFunction("x","y") {
//			public double value(Variable v) {
//				double x = v.get("x");
//				double y = v.get("y");
//				if(Math.sqrt((x-3.0)*(x-3.0)+(y-2.5)*(y-2.5))<0.5)
//					return 2.5;
//				else
//					return 0.0;
//			}			
//		};
//		NodeList nodes = mesh.getNodeList();
//		for(int j=1;j<=nodes.size();j++) {
//			a0.set(j,
//					fa0.value(
//							Variable.createFrom(fa0, nodes.at(j), j) ));
//		}
		plotVector(mesh, a0, "a0.dat");
		
		//v0: 求解状态方程得到
		Vector[] v0 = new SparseVector[nDataBlock];
		//lambda0: 求解伴随方程得到
		Vector[] lambda0 = new SparseVector[nDataBlock];
		
		for(int i=0;i<nDataBlock;i++) {
			ParamOfLightSource param =  paramList.get(i);
			v0[i] = this.solveStateEquation(param.lnu0_x,param.lnu0_y,a0,param.g_tidle);
			plotVector(mesh, v0[i], String.format("M%02d_v0.dat",i));
			lambda0[i] = this.solveAdjointEquation(param.s_i, param.lnu0_x, param.lnu0_y, 
				param.laplace_lnu0, a0, v0[i], param.g_tidle);
			plotVector(mesh, lambda0[i], String.format("M%02d_lambda0.dat",i));
		}
		//************************************************************
		
		//迭代求解
		BlockVector f = new SparseBlockVector(nDataBlock*2+1);
		for(int iter=0;iter<1;iter++) {
			
			this.iterNum = iter;
			
			BlockMatrix BM = this.getSearchDirectionMatrix(a0, v0, paramList);
			for(int i=1;i<=nDataBlock;i++) {
				ParamOfLightSource param =  paramList.get(i-1);
				//状态方程余量
				Vector resLlmd = this.getResLlmd(param.lnu0_x,  param.lnu0_y,  
						a0, v0[i-1],param.g_tidle).scale(-1.0);
		        plotVector(mesh,resLlmd,
		        		String.format("M%02d_Res_Llambda%02d.dat",i-1,this.iterNum));
				f.setBlock(nDataBlock+i, resLlmd);			
				
				//伴随方程余量
				Vector resLu = this.getResLu(param.s_i, 
					param.lnu0_x, param.lnu0_y, param.laplace_lnu0,
					a0, v0[i-1], param.g_tidle,lambda0[i-1]).scale(-1.0);//？1.0
		        plotVector(mesh,resLu,
		        		String.format("M%02d_Res_Lv%02d.dat",i-1,this.iterNum));

				f.setBlock(i, resLu);
			}
			//正则化参数方程余量
			f.setBlock(nDataBlock*2+1, this.getResLq(aGlob, 
					v0, lambda0, a0, 
					FC.c(aBkSquare) //a(x)的边界条件: back ground of a(x)
					).scale(-1.0));

			
			//设置边界条件并求解
			//需要交换矩阵BM的第1, 2列，然后设置边界条件，这样使得矩阵A, A'可逆
			BlockMatrix newBM = this.changeBlockColumn(BM, nDataBlock);
			this.imposeDirichletCondition(newBM, f, nDataBlock, FC.c0);
			SchurComplementLagrangianSolver solver = 
				new SchurComplementLagrangianSolver(BM, f,mesh);
			BlockVector x = solver.solveMulti();

			//this.imposeDirichletCondition(BM, f, FC.c0);
			//this.imposeDirichletCondition(BM, f, nDataBlock,FC.c0);
			//SolverJBLAS sol = new SolverJBLAS();
			//BlockVector x = (BlockVector)sol.solveDGESV(BM, f);
			
			for(int i=1;i<=nDataBlock;i++) {
				plotVector(mesh, x.getBlock(i), String.format("M%02d_delta_v%02d.dat",i-1,iter));
				plotVector(mesh, x.getBlock(nDataBlock+i), String.format("M%02d_delta_lambda%02d.dat",i-1,iter));
			}
			
			Vector delta_a = x.getBlock(nDataBlock*2+1);
				
			//!!!只截取delta_a还不够，还需要重新计算v和lambda
//			for(int j=1;j<=nodes.size();j++) {
//				Node node = nodes.at(j);
//				if(node.coord(1)<2 || node.coord(1)>4)
//					delta_a.set(j, 0.0);
//				else if(node.coord(2)<1.5)
//					delta_a.set(j, 0.0);
//					
//			}
			
			plotVector(mesh, delta_a, String.format("delta_a%02d.dat",iter));
			
			for(int i=0;i<nDataBlock;i++) {
				ParamOfLightSource param =  paramList.get(i);
				Vector vv2 = solveDeltaV(a0, param.lnu0_x, param.lnu0_y, 
						delta_a, v0[i]);
				plotVector(mesh, vv2, String.format("M%02d_rlt_v_2solve%02d.dat",i,iterNum));
			}
						
			
			//待定，与beta选取有关
			//double stepLength = 0.9*this.beta*10;
			//double stepLength = 0.0002;
			double stepLength = 0.9;
			
			for(int i=1;i<=nDataBlock;i++) {
				v0[i-1].add(stepLength, x.getBlock(i));
				lambda0[i-1].add(stepLength, x.getBlock(nDataBlock+i));
			}
			a0.add(stepLength, x.getBlock(nDataBlock*2+1));
			
			for(int i=1;i<=nDataBlock;i++) {
				plotVector(mesh, v0[i-1], String.format("M%02d_rlt_v%02d.dat",i-1,iter));
				plotVector(mesh, lambda0[i-1], String.format("M%02d_rlt_lambda%02d.dat",i-1,iter));
			}
			plotVector(mesh, a0, String.format("rlt_a%02d.dat",iter));
			
			//计算新a0对应的v
			for(int i=0;i<nDataBlock;i++) {
				ParamOfLightSource param =  paramList.get(i);
				Vector vv = solveStateEquation(param.lnu0_x,param.lnu0_y,a0,param.g_tidle);
				plotVector(mesh, vv, String.format("M%02d_rlt_v_solve%02d.dat",i,iterNum));
			}
		}
		return a0;

		
	}
	
	public void testRefine() {
		//ElementList eToRefine = computeRefineElement(mesh, alpha_avg_smooth, 0.06);
		ElementList eToRefine = new ElementList();
		eToRefine.add(mesh.getElementList().at(1));
		eToRefine.add(mesh.getElementList().at(2));
		eToRefine.add(mesh.getElementList().at(3));
		eToRefine.add(mesh.getElementList().at(15));
		eToRefine.add(mesh.getElementList().at(16));
		eToRefine.add(mesh.getElementList().at(17));
		this.refineMesh(eToRefine);
	}
	
	
	public void refine(Vector ak) {
		ElementList eToRefine = computeRefineElement(mesh, ak, 0.06);
		this.refineMesh(eToRefine);
	}	
	
	
	public static void main(String[] args) {
		VariationGaussNewtonDOT vgn = new VariationGaussNewtonDOT();
		vgn.debug = true;
//		vgn.readMeshTriangle();
		vgn.readMeshRectangle();
		
		
		//vgn.testRefine();
		for(int k=0;k<2;k++) {
			List<ParamOfLightSource> paramList = 
							new ArrayList<ParamOfLightSource>();
			NodeList nodes = vgn.mesh.getNodeList();
			
			for(int i=0; i<vgn.LS.length; i++) {
				//初始化模型“光源位置” 和 “包含物位置”
				vgn.reinitModelLight(i);
				
				ParamOfLightSource para = new ParamOfLightSource();
				para.s_i = i;
				para.u0 = vgn.solveU0(i);
				para.lnu0 = FMath.log(para.u0);
				para.lnu0_x = Tools.computeDerivative(vgn.mesh, para.lnu0, "x");
				para.lnu0_y = Tools.computeDerivative(vgn.mesh, para.lnu0, "y");
				para.laplace_lnu0 = Tools.computeLaplace2D(vgn.mesh, para.lnu0);
				
				para.lnu0_x = Utils.gaussSmooth(vgn.mesh, para.lnu0_x, 1, 0.5);
				para.lnu0_x = Utils.gaussSmooth(vgn.mesh, para.lnu0_x, 1, 0.5);
				para.lnu0_x = Utils.gaussSmooth(vgn.mesh, para.lnu0_x, 1, 0.5);
				para.lnu0_y = Utils.gaussSmooth(vgn.mesh, para.lnu0_y, 1, 0.5);
				para.lnu0_y = Utils.gaussSmooth(vgn.mesh, para.lnu0_y, 1, 0.5);
				para.lnu0_y = Utils.gaussSmooth(vgn.mesh, para.lnu0_y, 1, 0.5);
				para.laplace_lnu0 = Utils.gaussSmooth(vgn.mesh, para.laplace_lnu0, 1, 0.5);
				para.laplace_lnu0 = Utils.gaussSmooth(vgn.mesh, para.laplace_lnu0, 1, 0.5);
				para.laplace_lnu0 = Utils.gaussSmooth(vgn.mesh, para.laplace_lnu0, 1, 0.5);
				
				plotVector(vgn.mesh, para.lnu0_x, String.format("M%02d_lnu0_x.dat",i));
				plotVector(vgn.mesh, para.lnu0_y, String.format("M%02d_lnu0_y.dat",i));
				plotVector(vgn.mesh, para.laplace_lnu0, String.format("M%02d_laplace_lnu0.dat",i));
				
				Vector vReal = vgn.solveRealV(para.lnu0_x,para.lnu0_y);
		        plotVector(vgn.mesh,vReal,String.format("M%02d_vReal.dat",i));
	
				para.g_tidle = vReal.copy();
				//for(int j=1;j<=nodes.size();j++) {
				//	if(nodes.at(j).isInnerNode())
				//		para.g_tidle.set(j,0.0);
				//}
				plotVector(vgn.mesh, para.g_tidle, String.format("M%02d_g_tidle.dat",i));
				
				paramList.add(para);
			}
			
	
			//a(x)参考值（GCM方法得到的结果）
			vgn.aGlob = new SparseVector(nodes.size());
			for(int j=1;j<=nodes.size();j++) {
				double mu_a = vgn.modelGuess.mu_a.value(
								Variable.createFrom(vgn.modelGuess.mu_a, nodes.at(j), j));
				vgn.aGlob.set(j, mu_a/vgn.model_k);
			}
			VariationGaussNewtonDOT.plotVector(vgn.mesh, vgn.aGlob, "aGlob.dat");
	
			//vgn.gaussNewtonIterate(paramList);
			Vector ak = vgn.gaussNewtonIterateMulti(paramList);
			
			vgn.refine(ak);
		
		}
		
	}
	
	
	
	public void assignLinearShapFunction(Mesh mesh) {
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
	
	public ElementList computeRefineElement(Mesh mesh, Vector v, double persent) {
	    Vector v_smooth = Utils.gaussSmooth(mesh, v, 1, 0.5);

		ElementList eList = mesh.getElementList();
		ElementList eToRefine = new ElementList();

		List<PairDoubleInteger> list = new ArrayList<PairDoubleInteger>();
		for(int i=1;i<=eList.size();i++) {
			Element e = eList.at(i);
			Vector v1 = new SparseVector(e.nodes.size());
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
		return eToRefine;
	}	
}
