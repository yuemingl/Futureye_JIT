package edu.uta.futureye.application;

import java.io.File;
import java.util.ArrayList;
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
import edu.uta.futureye.core.NodeType;
import edu.uta.futureye.core.geometry.GeoEntity;
import edu.uta.futureye.function.AbstractFunction;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.DuDn;
import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.basic.Vector2Function;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.function.operator.FMath;
import edu.uta.futureye.io.MeshReader;
import edu.uta.futureye.io.MeshWriter;
import edu.uta.futureye.lib.assembler.AssemblerScalar;
import edu.uta.futureye.lib.element.FEBilinearRectangle;
import edu.uta.futureye.lib.element.FELinearTriangle;
import edu.uta.futureye.lib.weakform.WeakFormL22D;
import edu.uta.futureye.lib.weakform.WeakFormLaplace2D;
import edu.uta.futureye.util.FutureyeException;
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
 * v_tidle = ln{u}
 * 
 * @author liuyueming
 *
 */
public class VariationGaussNewtonDOT2 {
	public Mesh mesh;
	public Mesh meshBig;
	
	
	protected static String outputFolder = "Lagrangian_Klibanov2";
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
    double abk = 0.1/model_k;//back ground of a(x)
    
    //正则化参数
    //double beta = 0.01;
    //double beta = 1.0;
    //double beta = 10.0;
    double beta = 1000;

    int iterNum = 0;
    
    //光源x坐标位置数组
	double[] LS;

    /**
     * 构造函数，设定包含物位置
     */
    public VariationGaussNewtonDOT2() {
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
		//有包含物mu_a，真实模型
		modelReal.setMu_a(2.91, 2.36, 0.5,
				2.0, //peak value of mu_a
				1); //Number of inclusions
		//有包含物mu_a，猜测模型
		modelGuess.setMu_a(2.91, 2.36, 0.5,
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
	
	/**
	 * Solve background forward problem u0
	 * 
	 * @param s_i
	 * @return
	 */
	public Vector solveU0(int s_i) {
		//在mesh上求解，top边界条件不好确定，因此在比较大的meshBig上求解，
		//然后截取出来到较小的mesh上
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
	
	public Vector solveGuessU(int s_i) {
		Vector u0Big = modelGuess.solveNeumann(meshBig);
		plotVector(meshBig, u0Big, "u0_big"+s_i+".dat");
		Vector u0Samll = Tools.extractData(meshBig, mesh, u0Big);
		plotVector(mesh, u0Samll, "u0_"+s_i+".dat");
		return u0Samll;
	}	
	
	/**
	 * 
	 * @param s_i
	 * @param g_tidle = realU|_\Gamma
	 * @return
	 */
	public Vector solveGuessU(int s_i, Vector g_tidle) {
		Vector u0 = modelGuess.solveDirichlet(mesh, new Vector2Function(g_tidle));
		plotVector(mesh, u0, "u0_"+s_i+".dat");
		return u0;
	}	
	
	public Vector solveRealU(int s_i) {
		Vector u0Big = modelReal.solveNeumann(meshBig);
		plotVector(meshBig, u0Big, "u0_big"+s_i+".dat");
		Vector u0Samll = Tools.extractData(meshBig, mesh, u0Big);
		plotVector(mesh, u0Samll, "u0_"+s_i+".dat");
		return u0Samll;
	}	
	
	/**
	 * Get stiff matrix and load vector for equation of lambda
	 * 
	 * (v-z, dv) - (\nabla{dv},\nabla{lmd}) + (2*\nabla{v}\cdot\nabla{dv},lmd)=0
	 * 
	 * lmd=0 in \partial\Omega
	 * 
	 * @param s_i
	 * @param v
	 * @return
	 */
	public Equation getEqnLambda(int s_i,  
			Vector v, Vector z) {
		Vector v_x = Tools.computeDerivative(mesh, v, "x");
		Vector v_y = Tools.computeDerivative(mesh, v, "x");
		Function b1 = new Vector2Function(v_x);
		Function b2 = new Vector2Function(v_y);
		WeakFormGCMDual weakForm = new WeakFormGCMDual();

		Vector v_g = FMath.axpy(-1.0, z, v);
		Function fv_g = new Vector2Function(v_g);
		plotFunction(mesh, fv_g, String.format("v_g%02d.dat",s_i));

		weakForm.setF(fv_g.M(FC.c(-1.0)));

		weakForm.setParam(
				FC.c(-1.0),//注意，有负号!!!
				FC.c0, 
				FC.c(2.0).M(b1),
				FC.c(2.0).M(b2)
			);
		
		mesh.clearBorderNodeMark();
		HashMap<NodeType, Function> mapNTF = new HashMap<NodeType, Function>();
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
		eqn.f = load;
		
		return eqn;
	}	

	
	public static void main(String[] args) {
		VariationGaussNewtonDOT2 vgn = new VariationGaussNewtonDOT2();
		vgn.debug = true;
		vgn.readMeshTriangle();
//		vgn.readMeshRectangle();

		List<ParamOfLightSource> paramList = 
						new ArrayList<ParamOfLightSource>();
		
		for(int i=0; i<vgn.LS.length; i++) {
			//初始化模型“光源位置” 和 “包含物位置”
			vgn.reinitModelLight(i);
			
			//ParamOfLightSource para = new ParamOfLightSource();
			//para.s_i = i;
			//para.u0 = vgn.solveU0(i);
			//para.lnu0 = FMath.log(para.u0);

			Vector realU = vgn.solveRealU(i);
			Vector guessU = vgn.solveGuessU(i,realU);
	        plotVector(vgn.mesh,guessU,String.format("M%02d_guessU.dat",i));
	        plotVector(vgn.mesh,realU,String.format("M%02d_realU.dat",i));

	        Vector z = FMath.log(realU);
	        Vector v = FMath.log(guessU);
	        plotVector(vgn.mesh,z,String.format("M%02d_z.dat",i));
	        plotVector(vgn.mesh,v,String.format("M%02d_v.dat",i));
	        
	        Equation eq = vgn.getEqnLambda(i, v, z);
	        Solver sol = new Solver();
	        Vector lmd = sol.solveCGS(eq.A, eq.f);
	        plotVector(vgn.mesh,lmd,String.format("M%02d_lambda.dat",i));
	        
	        
			//paramList.add(para);
		}
	}
}
