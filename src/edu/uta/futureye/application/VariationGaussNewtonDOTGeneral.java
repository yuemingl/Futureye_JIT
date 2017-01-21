package edu.uta.futureye.application;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import edu.uta.futureye.algebra.SparseBlockMatrix;
import edu.uta.futureye.algebra.SparseBlockVector;
import edu.uta.futureye.algebra.SparseMatrixRowMajor;
import edu.uta.futureye.algebra.SparseVectorHashMap;
import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.SparseMatrix;
import edu.uta.futureye.algebra.intf.SparseVector;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.algebra.solver.external.SolverJBLAS;
import edu.uta.futureye.core.DOF;
import edu.uta.futureye.core.DOFOrder;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.core.NodeRefined;
import edu.uta.futureye.core.NodeType;
import edu.uta.futureye.core.Refiner;
import edu.uta.futureye.core.geometry.GeoEntity;
import edu.uta.futureye.function.MultiVarFunc;
import edu.uta.futureye.function.FMath;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.DuDn;
import edu.uta.futureye.function.basic.DuDx;
import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.basic.Vector2MathFunc;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.io.MeshReader;
import edu.uta.futureye.io.MeshWriter;
import edu.uta.futureye.lib.assembler.AssemblerScalar;
import edu.uta.futureye.lib.element.FEBilinearRectangleRegular;
import edu.uta.futureye.lib.element.FELinearTriangleOld;
import edu.uta.futureye.lib.weakform.WeakFormLaplace2D;
import edu.uta.futureye.util.FutureyeException;
import edu.uta.futureye.util.PropertiesReader;
import edu.uta.futureye.util.Utils;
import edu.uta.futureye.util.container.DOFList;
import edu.uta.futureye.util.container.ElementList;
import edu.uta.futureye.util.container.NodeList;
import edu.uta.futureye.util.container.ObjIndex;

/**
 * Implementation of the paper: 
 *   'A Framework For The Adaptive Finite Element Solution Of Large-Scale Inverse Problems'
 * 
 * Lagrange Multiplier Method based on the model:
 * 
 * -\nabla{(1/(a*k))*\nabla{u}} + u = \delta/a
 * 
 * where 
 *   k = 3*mu_s'
 *   a = a(x) = mu_a(x)
 * 
 * Measurements take place in the whole domain of \Omega or on the boundary of \Gamma,
 * where \Gamma = \partial\Omega 
 *
 * @author liuyueming
 *
 */
public class VariationGaussNewtonDOTGeneral {
	public boolean debug = false;
	
	//输出目录（格式：outputFolderBase+outputFolderIndex）
	protected String outputFolderBase;
	protected int outputFolderIndex = 0;
	
	//网格
	public Mesh mesh;
	public Mesh meshBig;
	public String gridFileBig;
	public String gridFileSmall;
	
	//是否使用向量 mu_a，由于向量与网格有关，因此在不同网格上求解问题需要对mu_a插值 
	boolean useVectorMu_a = false;
	//已知参数（向量类型的mu_a）
	protected Vector aReal = null;
	protected Vector aGuess = null;
	protected Vector aInit = null;
	protected Mesh aMesh = null;
	public String inputDataFolder;
	public String inputDataFile_aReal;
	public String inputDataFile_aGuess;
	public String inputDataMesh;

	//k = 3*mu_s'
	MathFunc model_k = FC.c(50.0);
	//Background of a(x) = mu_a(x) = 0.1
	double aBackground = 0.1;

	//模型定义
	ModelDOTMult modelReal = new ModelDOTMult(); //Real inclusion model
	ModelDOTMult modelGuess = new ModelDOTMult();//Guess model from GCM of inclusion 
	ModelDOTMult modelInit = new ModelDOTMult(); //Initial model of inclusion

	//测量类型开关，以下为u_g在整个区域上都已知
	public boolean bTestWholdDomain = true; //=true测量数据在整个区域，=false测量数据在边界上 （=true 有强烈震荡，为什么？需要调整beta较大些可以解决！）
	public boolean bTestWholeDomainDirichletBoundary = true; //测量数据在整个区域并且=ture时使用Dirichlet边界条件，=false时使用Neumann边界条件
	public boolean bTestBoundaryAsWholdDomain = false;
    
    //mu_a from GCM
	protected Vector aGlob = null;//GCM方法得到的a_glob(x)
	
    //正则化参数
    //double beta = 1.0; //bTestWholdDomain = true;
    //protected double beta = 0.03; //bTestWholdDomain = true;
    protected double beta = 0.006; //bTestWholdDomain = true;
    //整个区域
    //double beta = 100000; //bTestWholdDomain = false;
    
    //Current iterate number
    protected int iterNum = 0;
    //Current refine num
    protected int refineNum = 0;
    
    //Parameters of refinement and iteration
    public int totalRefineNum = 3;
    public double[] refineFactorsForMesh = null;
    public int[] maxIterNumPerRefinement = null;
    
    public double[] initStepLengthPerRefinement = null;
	public int maxSearchNum = 30; //最大搜索次数
	
	//利用以下因子控制计算结果
	public double stepReduceFactor = 0.75; //搜索时步长减小因子
	public double maxTargetNormIncFactor = 3.0; //迭代一步后，目标范数增加因子
	public double maxInfNormIncFactor = 1.2; //迭代一步后，最大范数增加因子
	
	//For every delta_a
	public double[] cutX = {-100,100}; //big enough = no cut in X direction
	public double[] cutY = {-100,100}; //big enough = no cut in Y direction
	public double cutThreshold = 0.5; //[0.0,1.0]
	public int smoothNum = 4;// >=0
	
    //光源\vec{x}坐标位置数组
    protected double[] LSx = null;
    protected double[] LSy = null;

    protected FileOutputStream out = null;
    protected PrintWriter br = null;
	
    //对应每个光源s_i的参数，包括测量数据
    public static class ParamOfLightSource {
    	int s_i;
    	Vector g; //= u|_\Gamma
    }
    
	/**
	 * Initialization of configuration parameters
	 * 
	 */
	public void init() {
		
//test9		
//		//有包含物mu_a，真实模型
//		modelReal.setMu_a(2.40, 2.60, 0.6,
//				0.4, //peak value of mu_a
//				1); //Number of inclusions
//		//有包含物mu_a，猜测模型
//		modelGuess.setMu_a(2.40, 2.60, 0.6,
//				0.2, //peak value of mu_a
//				1); //Number of inclusions
//		//有包含物mu_a，迭代初始值
//		modelInit.setMu_a(2.40, 2.60, 0.6,
//				0.2, //peak value of mu_a
//				1); //Number of inclusions

//test9_1
//		//有包含物mu_a，真实模型
//		modelReal.setMu_a(2.20, 2.60, 0.6,
//				0.4, //peak value of mu_a
//				1); //Number of inclusions
//		//有包含物mu_a，猜测模型
//		modelGuess.setMu_a(2.40, 2.60, 0.6,
//				0.4, //peak value of mu_a
//				1); //Number of inclusions
//		//有包含物mu_a，迭代初始值
//		modelInit.setMu_a(2.40, 2.60, 0.6,//0.8
//				0.4, //peak value of mu_a
//				1); //Number of inclusions

//test9_2
//		//有包含物mu_a，真实模型
//		modelReal.setMu_a(2.40, 2.40, 0.6,
//				0.4, //peak value of mu_a
//				1); //Number of inclusions
//		//有包含物mu_a，猜测模型
//		modelGuess.setMu_a(2.40, 2.60, 0.6,
//				0.4, //peak value of mu_a
//				1); //Number of inclusions		

//test9_3
//		//有包含物mu_a，真实模型
//		modelReal.setMu_a(2.40, 2.60, 0.6,
//				0.4, //peak value of mu_a
//				1); //Number of inclusions
//		//有包含物mu_a，猜测模型
//		modelGuess.setMu_a(2.40, 2.40, 0.6,
//				0.4, //peak value of mu_a
//				1); //Number of inclusions		
		
		
//		//有包含物mu_a，真实模型
//		modelReal.setMu_a(2.20, 2.0, 0.6,
//				0.4, //peak value of mu_a
//				1); //Number of inclusions
//		//有包含物mu_a，猜测模型
//		modelGuess.setMu_a(2.40, 2.0, 0.6,
//				0.4, //peak value of mu_a
//				1); //Number of inclusions
//		//有包含物mu_a，迭代初始值
//		modelInit.setMu_a(2.40, 2.0, 0.6,//0.8
//				0.4, //peak value of mu_a
//				1); //Number of inclusions
		
		//有包含物mu_a，真实模型
		//2011-10-11
//		modelReal.setMu_a(2.50, 2.0, 0.6,
//				0.4, //peak value of mu_a
//				1); //Number of inclusions
//		
//		//有包含物mu_a，猜测模型
//		modelGuess.setMu_a(2.45, 2.0, 0.6,
//				0.36, //peak value of mu_a
//				1); //Number of inclusions
//		//有包含物mu_a，迭代初始值
//		modelInit.setMu_a(2.45, 2.0, 0.6,//0.8
//				0.36, //peak value of mu_a
//				1); //Number of inclusions

		//自己构造的mu_a
//		modelReal.mu_a = generateTestRealMu_a(0.4, 0.1);
//		modelGuess.mu_a = generateTestGuessMu_a(0.4,0.1);
//		modelInit.mu_a = generateTestGuessMu_a(0.4,0.1);
		
		
		
//		//自己构造的mu_a（通过class ModelPoissonEx.gen1()），更光滑，在整个区域都有值
//		useVectorMu_a = true;
//		aMesh = mesh.copy();
//		aReal = DataReader.readVector(String.format("./"+this.getOutputFolder()+"/Input/input_real_mu_a.dat"));
//		aGuess = DataReader.readVector(String.format("./"+this.getOutputFolder()+"/Input/input_guess_mu_a.dat"));
//		aInit = aGuess.copy();
		
//		//test
//		//Function faReal = generateTestRealMu_a(0.4, 0.1);
//		//Function faGuess = generateTestGuessMu_a(0.4,0.1);
//		//aReal = Tools.function2vector(mesh, faReal);
//		//aGuess = Tools.function2vector(mesh, faGuess);		
//		//aInit = aGuess.copy();
		
		//自己构造的mu_a（通过class ModelPoissonEx.gen2()），两个inclusion
//		useVectorMu_a = true;
//		aMesh = mesh.copy();
//		aReal = DataReader.readVector(String.format("./"+this.getOutputFolder()+"/Input/input_real_mu_a_test2v2.dat"));
//		aGuess = DataReader.readVector(String.format("./"+this.getOutputFolder()+"/Input/input_guess_mu_a_test2.dat"));
//		aInit = aGuess.copy();
		
//		//test15	
		//有包含物mu_a，真实模型
//		modelReal.setMu_a(ModelDOTMult.getMu_a(0.05, 0.0, 0.2,
//				0.3, //peak value of mu_a
//				1)); //Number of inclusions
//		//有包含物mu_a，猜测模型
//		modelGuess.setMu_a(ModelDOTMult.getMu_a(0.0, 0.0, 0.2,
//				0.3, //peak value of mu_a
//				1)); //Number of inclusions
//		//有包含物mu_a，迭代初始值
//		modelInit.setMu_a(ModelDOTMult.getMu_a(0.0, 0.0, 0.2,
//				0.3, //peak value of mu_a
//				1)); //Number of inclusions
//		plotFunction(mesh,modelGuess.getMu_a(),"mu_a_function_guess.dat");
//		plotFunction(mesh,modelReal.getMu_a(),"mu_a_function_real.dat");
		
		
		useVectorMu_a = true;
		aMesh = mesh.copy();
		aReal = DataReader.readVector(String.format(inputDataFolder+inputDataFile_aReal));
		aGuess = DataReader.readVector(String.format(inputDataFolder+inputDataFile_aGuess));
		//如果网格相同，插值这一步可以省略
		if(inputDataMesh != null && !inputDataMesh.isEmpty()) {
			MeshReader mReader = new MeshReader(inputDataFolder+inputDataMesh);
			Mesh inputMesh = mReader.read2DMesh();
			aReal = Tools.interplateFrom(inputMesh, mesh, new Vector2MathFunc(aReal,inputMesh,"x","y"));
			aGuess = Tools.interplateFrom(inputMesh, mesh, new Vector2MathFunc(aGuess,inputMesh,"x","y"));
		}
        aInit = aGuess.copy();
	}
	
	public void readParameters(String configFileName) {
		PropertiesReader pReader = new PropertiesReader(configFileName);
		Boolean debug = pReader.getBoolean("debug");
		if(debug != null) this.debug = debug;

		//结果输出目录前缀（每次加密网格后创建一个新目录，结尾数字编号自动增加）
		String outputFolder = pReader.getString("outputFolder");
		if(outputFolder == null || outputFolder.isEmpty()) throw new FutureyeException("Please specify 'outputFolder' parameter in config file!");
		this.outputFolderBase = outputFolder;
		String oFolder = getOutputFolder();
	    if(!oFolder.isEmpty()) {
		    File file = new File(oFolder);
			if(!file.exists()) {
				file.mkdirs();
			}
	    }
		
		String gridFileBig = pReader.getString("gridFileBig");
		if(gridFileBig == null || gridFileBig.isEmpty()) throw new FutureyeException("Please specify 'gridFileBig' parameter in config file!");
		this.gridFileBig = gridFileBig;
		String gridFileSmall = pReader.getString("gridFileSmall");
		if(gridFileSmall == null || gridFileSmall.isEmpty()) throw new FutureyeException("Please specify 'gridFileSmall' parameter in config file!");
		this.gridFileSmall = gridFileSmall;    	
		
		//Light sources (positions)
		double[] LSx = pReader.getDoubleArray("LSx");
		if(LSx != null) this.LSx = LSx;
		double[] LSy = pReader.getDoubleArray("LSy");
		if(LSy != null) this.LSy = LSy;
		
		//Mesh refinement control
		Integer totalRefineNum = pReader.getInteger("totalRefineNum");
		if(totalRefineNum != null) this.totalRefineNum = totalRefineNum;
		double[] refineFactorsForMesh = pReader.getDoubleArray("refineFactorsForMesh");
		if(refineFactorsForMesh != null) this.refineFactorsForMesh = refineFactorsForMesh;
		int[] maxIterNumPerRefinement = pReader.getIntegerArray("maxIterNumPerRefinement");
		if(maxIterNumPerRefinement != null) this.maxIterNumPerRefinement = maxIterNumPerRefinement;
		double[] initStepLengthPerRefinement = pReader.getDoubleArray("initStepLengthPerRefinement");
		if(initStepLengthPerRefinement != null) this.initStepLengthPerRefinement = initStepLengthPerRefinement;
		
		//Iteration control
		Double stepReduceFactor = pReader.getDouble("stepReduceFactor");
		if(stepReduceFactor != null) this.stepReduceFactor = stepReduceFactor;
		Double maxTargetNormIncFactor = pReader.getDouble("maxTargetNormIncFactor");
		if(maxTargetNormIncFactor != null) this.maxTargetNormIncFactor = maxTargetNormIncFactor;
		Double maxInfNormIncFactor = pReader.getDouble("maxInfNormIncFactor");
		if(maxInfNormIncFactor != null) this.maxInfNormIncFactor = maxInfNormIncFactor;
		
		//Cut control for every delta_a
		double[] cutX = pReader.getDoubleArray("cutX");
		if(cutX != null) this.cutX = cutX;
		double[] cutY = pReader.getDoubleArray("cutY");
		if(cutY != null) this.cutY = cutY;
		Double cutThreshold = pReader.getDouble("cutThreshold");
		if(cutThreshold != null) this.cutThreshold = cutThreshold;
		Integer smoothNum = pReader.getInteger("smoothNum");
		if(smoothNum != null) this.smoothNum = smoothNum;
		
		//输入数据文件aReal aGuess aMesh
		String inputDataFolder = pReader.getString("inputDataFolder");
		if(inputDataFolder == null || inputDataFolder.isEmpty()) throw new FutureyeException("Please specify 'inputDataFolder' parameter in config file!");
		this.inputDataFolder = inputDataFolder;
		String inputDataFile_aReal = pReader.getString("inputDataFile_aReal");
		this.inputDataFile_aReal = inputDataFile_aReal;
		String inputDataFile_aGuess = pReader.getString("inputDataFile_aGuess");
		if(inputDataFile_aGuess == null || inputDataFile_aGuess.isEmpty()) throw new FutureyeException("Please specify 'inputDataFile_aGuess' parameter in config file!");
		this.inputDataFile_aGuess = inputDataFile_aGuess;
		String inputDataMesh = pReader.getString("inputDataMesh");
		this.inputDataMesh = inputDataMesh;

		//正则化参数
		Double beta = pReader.getDouble("beta");
		if(beta != null) this.beta = beta;

	}
	
//	/**
//	 * Override some parameters by args
//	 * @param args
//	 */
//	public void overrideParameters(String[] args) {
//		if(args.length == 3) {
//			this.stepReduceFactor = Double.parseDouble(args[0]);
//			this.maxTargetNormIncFactor = Double.parseDouble(args[1]);
//			this.maxInfNormIncFactor = Double.parseDouble(args[2]);
//			System.out.println("---begin with args specified, override configure file!---");
//		}
//	}
	
	public String getString(double[] ary) {
		if(ary == null || ary.length == 0) return "";
		StringBuilder sb = new StringBuilder();
		sb.append(ary[0]);
		for(int i=1;i<ary.length;i++) {
			sb.append(",");
			sb.append(ary[i]);
		}
		return sb.toString();
	}
	public String getString(int[] ary) {
		if(ary == null || ary.length == 0) return "";
		StringBuilder sb = new StringBuilder();
		sb.append(ary[0]);
		for(int i=1;i<ary.length;i++) {
			sb.append(",");
			sb.append(ary[i]);
		}
		return sb.toString();
	}	
	
	public void printParameters() {
		System.out.println(String.format("---------parameters----------\n"+
				"debug=%s\n"+
				"\n"+
				
				"outputFolder=%s\n"+
				"gridFileBig=%s\n"+
				"gridFileSmall=%s\n"+
				"\n"+
				
				"LSx=%s\n"+
				"LSy=%s\n"+
				"\n"+
				
				"totalRefineNum=%s\n"+
				"refineFactorsForMesh=%s\n"+
				"maxIterNumPerRefinement=%s\n"+
				"initStepLengthPerRefinement=%s\n"+
				"\n"+
				
				"stepReduceFactor=%s\n"+
				"maxTargetNormIncFactor=%s\n"+
				"maxInfNormIncFactor=%s\n"+
				"\n"+
				
				"cutX=%s\n"+
				"cutY=%s\n"+
				"cutThreshold=%s\n"+
				"smoothNum=%s\n"+
				"\n"+
				
				"inputDataFolder=%s\n"+
				"inputDataFile_aReal=%s\n"+
				"inputDataFile_aGuess=%s\n"+
				"inputDataMesh=%s\n"+
				"\n"+
				
				"beta=%s\n",

				this.debug,
				this.outputFolderBase,
				this.gridFileBig,
				this.gridFileSmall,

				getString(this.LSx),
				getString(this.LSy),

				this.totalRefineNum,
				getString(this.refineFactorsForMesh),
				getString(this.maxIterNumPerRefinement),
				getString(this.initStepLengthPerRefinement),
				
				this.stepReduceFactor,
				this.maxTargetNormIncFactor,
				this.maxInfNormIncFactor,
				
				getString(this.cutX),
				getString(this.cutY),
				this.cutThreshold,
				this.smoothNum,
				
				this.inputDataFolder,
				this.inputDataFile_aReal,
				this.inputDataFile_aGuess,
				this.inputDataMesh,
				
				this.beta
				));
		System.out.println("-------------------------------\n");		
	}
	
    public String getOutputFolder() {
    	return String.format(outputFolderBase+"%02d", outputFolderIndex);
    }
    
	public void setOutputFolderIndex(int index) {
		outputFolderIndex = index;
	}	
	
	public static MathFunc generateTestRealMu_a(double max, double bk) {
		final double fmax = max;
		final double fbk = bk;
		return new MultiVarFunc("x","y"){
			@Override
			public double apply(Variable v) {
				double x = v.get("x");
				double y = v.get("y");
				if(x>=1.9 && x<=2.2 && y>=2.1 && y<=2.52)
					return fmax+Math.random()*0.01;
				else if(x>=1.75 && x<=2.2 && y>=1.6 && y<=2.1)
					return fmax+Math.random()*0.01;
				else if(x>=1.9 && x<=2.2 && y>=1.45 && y<=1.6)
					return fmax+Math.random()*0.01;
				else
					return fbk;
			}
		};
	}
	
	public static MathFunc generateTestGuessMu_a(double max, double bk) {
		final double fmax = max;
		final double fbk = bk;
		return new MultiVarFunc("x","y"){
			@Override
			public double apply(Variable v) {
				double x = v.get("x");
				double y = v.get("y");
				if(x>=1.9 && x<=2.2 && y>=1.4 && y<=2.52)
					return fmax+Math.random()*0.01;
				//else if(x>=2.2 && x<=2.4 && y>=1.6 && y<=2.4)
				//	return fmax+Math.random()*0.01;
				else
					return fbk;
			}
		};
	}
	
	public static Vector2MathFunc generateTestRealMu_a2(Mesh mesh,double bk) {
		int[] nodes = {134,135,136,151,152,153,168,169,170,185,186,187};
		double[] values = {0.25,0.3,0.22,0.23,0.4,0.2,0.3,0.4,0.22,0.23,0.25,0.24};
		int[] nodes2 = {137,138,154,155,171,172};
		double[] values2 = {0.21,0.24,0.4,0.35,0.4,0.22};
		NodeList nodeList = mesh.getNodeList();
		int dim = nodeList.size();
		SparseVector v = new SparseVectorHashMap(dim);
		for(int i=1;i<=nodeList.size();i++) {
			v.set(i,bk);
			for(int j=0;j<nodes.length;j++) {
				if(i==nodes[j])
					v.set(i, values[j]);
			}
			for(int j=0;j<nodes2.length;j++) {
				if(i==nodes2[j])
					v.set(i, values2[j]);
			}
		}
		return new Vector2MathFunc(v);
	}
	
	public static Vector2MathFunc generateTestGuessMu_a2(Mesh mesh,double bk) {
		int[] nodes = {135,136,137,151,152,153,154,168,169,170,171,186,187};
		double[] values = {0.32,0.3,0.32,0.3,0.35,0.28,0.38,0.25,0.3,0.25,0.3,0.3,0.3};

		NodeList nodeList = mesh.getNodeList();
		int dim = nodeList.size();
		SparseVector v = new SparseVectorHashMap(dim);
		for(int i=1;i<=nodeList.size();i++) {
			v.set(i,bk);
			for(int j=0;j<nodes.length;j++) {
				if(i==nodes[j])
					v.set(i, values[j]);
			}
		}
		return new Vector2MathFunc(v);
	}
	
	public static Vector2MathFunc generateRealMu_aTest15(Mesh mesh,double bk) {
		//v1
//		int[] nodes = {
//				113,
//			125,126,
//			138,139,
//			151,152};
//		double[] values = {
//				0.4,
//			0.4,0.4,
//			0.4,0.4,
//			0.3,0.3};
//		int[] nodes2 = {128,141,142,154};
//		double[] values2 = {0.45,0.45,0.3,0.4};
		//v2=guess只有一个  v3=有两个guess，尽量接近
//		int[] nodes2 = {
//				115,116,
//				128,129,130,
//				141,142,143,
//				154,155,156,
//				167,168};
//		double[] values2 = {
//				0.2,  0.3,
//				0.2,  0.45, 0.3,
//				0.2,  0.45, 0.3,
//				0.2,  0.4 , 0.3,
//				0.3,  0.3};
		
//		//v4=很小的两个
//		int[] nodes = {139};
//		double[] values = {0.4};
//		int[] nodes2 = {182};
//		double[] values2 = {bk};
		
		//v5
		int[] nodes = {140,141,153,154};
		double[] values = {0.4,0.4,0.4,0.4};
		int[] nodes2 = {};
		double[] values2 = {};
		
		NodeList nodeList = mesh.getNodeList();
		int dim = nodeList.size();
		SparseVector v = new SparseVectorHashMap(dim);
		for(int i=1;i<=nodeList.size();i++) {
			v.set(i,bk);
			for(int j=0;j<nodes.length;j++) {
				if(i==nodes[j])
					v.set(i, values[j]);
			}
			for(int j=0;j<nodes2.length;j++) {
				if(i==nodes2[j])
					v.set(i, values2[j]);
			}
		}
		return new Vector2MathFunc(v);
	}
	
	public static Vector2MathFunc generateGuessMu_aTest15(Mesh mesh,double bk) {
		//v1
//		int[] nodes = {
//				113,
//				126,127,128,
//			138,139,140,141,
//			    152,153,154};
//		double[] values = {
//				0.25,
//				0.4,0.3,0.25,
//			0.4,0.3,0.4,0.4,
//				0.3,0.3,0.3};
		//v2,v3
//		int[] nodes = {
//				113,114,
//			125,126,127,
//			138,139,140,
//			151,152,153};
//		double[] values = {
//				0.4,0.1,
//			0.4,0.4,0.4,
//			0.4,0.4,0.4,
//			0.3,0.3,0.1};
//		int[] nodes2 = {
//				115,116,
//				128,129,130,
//				141,142,143,
//				154,155,156,
//				167,168};
//		double[] values2 = {
//				0.2,  0.3,
//				0.2,  0.45, 0.3,
//				0.2,  0.45, 0.3,
//				0.2,  0.4 , 0.3,
//				0.3,  0.3};
		
//		//v4
//		int[] nodes = {140};
//		double[] values = {0.4};
//		int[] nodes2 = {182};
//		double[] values2 = {bk};
		
		//v5
		int[] nodes = {141,142,154,155};
		double[] values = {0.35,0.35,0.35,0.35};
		int[] nodes2 = {};
		double[] values2 = {};
		
		NodeList nodeList = mesh.getNodeList();
		int dim = nodeList.size();
		SparseVector v = new SparseVectorHashMap(dim);
		for(int i=1;i<=nodeList.size();i++) {
			v.set(i,bk);
			for(int j=0;j<nodes.length;j++) {
				if(i==nodes[j])
					v.set(i, values[j]);
			}
			for(int j=0;j<nodes2.length;j++) {
				if(i==nodes2[j])
					v.set(i, values2[j]);
			}
		}
		return new Vector2MathFunc(v);
	}
	
    /**
     * 初始化模型“光源位置”
     * @param s_i：Light source No. (0,1,2...)
     */
    public void reinitModelLight(int s_i) {
		modelReal.setLightPosition(LSx[s_i], LSy[s_i]);
		modelGuess.setLightPosition(LSx[s_i], LSy[s_i]);
		modelInit.setLightPosition(LSx[s_i], LSy[s_i]);
    }
    

	public void plotVector(Mesh mesh, Vector v, String fileName) {
		String folder = getOutputFolder();
		Tools.plotVector(mesh, folder, fileName, v);
	}

	public void plotFunction(Mesh mesh, MathFunc fun, String fileName) {
		String folder = getOutputFolder();
		Tools.plotFunction(mesh, folder, fileName, fun);
	}	
	

	/**
	 * Log information in the file output_log.txt in every different folder
	 * 
	 */
	public void beginLog() {
		String folder = getOutputFolder();
		try {
			File file = new File("./"+folder+"/output_log.txt");
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
	
	/**
	 * Read mesh and assign DOF to elements
	 * 
	 * Supports:
	 * 	triangle elements
	 *  rectangle elements
	 *  
	 */
	public void readMesh(String folder){
        MeshReader readerBig = new MeshReader(folder+"/"+gridFileBig);
        MeshReader readerSmall = new MeshReader(folder+"/"+gridFileSmall);
        meshBig = readerBig.read2DMesh();
        mesh = readerSmall.read2DMesh();
        meshBig.computeNodeBelongsToElements();
        mesh.computeNodeBelongsToElements();
        mesh.computeNeighborNodes();
      
        //2.Mark border types
        HashMap<NodeType, MathFunc> mapNTF =
                new HashMap<NodeType, MathFunc>();
        mapNTF.put(NodeType.Dirichlet, null);
        mesh.markBorderNode(mapNTF);
        
        
        //3.Use element library to assign degrees of freedom (DOF) to element
        ElementList eList = mesh.getElementList();
        FELinearTriangleOld feLT = new FELinearTriangleOld();
//        FEBilinearRectangle feBR = new FEBilinearRectangle();
        FEBilinearRectangleRegular feBRR = new FEBilinearRectangleRegular();
        for(int i=1;i<=eList.size();i++) {
        	Element e = eList.at(i);
        	if(e.nodes.size()%3 == 0)
        		feLT.assignTo(eList.at(i));
        	else if(e.nodes.size()%4 == 0)
        		feBRR.assignTo(eList.at(i));
        }
  
        ElementList eListBig = meshBig.getElementList();
		for(int i=1;i<=eListBig.size();i++) {
	        Element e = eListBig.at(i);
	    	if(e.nodes.size()%3 == 0)
	    		feLT.assignTo(eListBig.at(i));
	    	else if(e.nodes.size()%4 == 0)
	    		feBRR.assignTo(eListBig.at(i));
		}
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
		
		System.out.println("Before refine [meshBig]: Element="+meshBig.getElementList().size()+", Node="+meshBig.getNodeList().size());
		Refiner.refineOnce(meshBig, eToRefineBig);
		System.out.println("After  refine [meshBig]: Element="+meshBig.getElementList().size()+", Node="+meshBig.getNodeList().size());

		System.out.println("Before refine [mesh]   : Element="+mesh.getElementList().size()+", Node="+mesh.getNodeList().size());
		Refiner.refineOnce(mesh, eToRefine);
		System.out.println("After  refine [mesh]   : Element="+mesh.getElementList().size()+", Node="+mesh.getNodeList().size());
		
		Tools.assignLinearShapFunction(meshBig);
		Tools.assignLinearShapFunction(mesh);
	}
	
	public void constrainHangingNodes(Mesh mesh, Vector v) {
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
	
	public void zeroHangingNode(Mesh mesh, Vector v) {
		NodeList nodes = mesh.getNodeList();
		for(int i=1;i<=nodes.size();i++) {
			Node node = nodes.at(i);
			if(node instanceof NodeRefined) {
				NodeRefined nRefined = (NodeRefined)node;
				if(nRefined.isHangingNode()) {
					v.set(node.globalIndex, 0.0);
				}
			}
		}
	}
	
	/**
	 * Get stiff matrix and load vector for equation of 'lambda'
	 * 
	 * L_u(\lambda)=0:
	 *   ((1/(a*k))*\nabla{\psi},\nabla{\lambda}) + (\psi,\lambda) = -(u-g,\psi)_\Omega
	 * 
	 * Boundary Condition: Dirichlet
	 *   \lambda = 0 on \Omega
	 * 
	 * where 
	 *    \Gamma = \partial\Omega
	 *    g = uReal|_\Gamma
	 * 
	 * Note:
	 *    bTestWholdDomain == true
	 *      bTestWholeDomainDirichletBoundary == true:  测量整个区域，Dirichlet边界条件
	 *      bTestWholeDomainDirichletBoundary == false: 测量跟个区域，Neumann边界条件
	 *    bTestWholdDomain == false: 测量边界区域，Neumann边界条件
	 * 
	 * @param s_i
	 * @param u: uk
	 * @param g: uReal|_\Gamma
	 * @return
	 */
	public Equation getEqnLambda(int s_i, Vector a, Vector u, Vector g) {
		WeakFormLaplace2D weakForm = new WeakFormLaplace2D();
		//(u - g)_\Gamma
		Vector u_g = FMath.axpy(-1.0, g, u);
		
		NodeList nodes = mesh.getNodeList();
		if(!bTestWholdDomain) {
			for(int j=1;j<=nodes.size();j++) {
				if(nodes.at(j).isInnerNode())
					u_g.set(j,0.0);
			}
		}
		
		MathFunc fu_g = new Vector2MathFunc(u_g);
		plotFunction(mesh, fu_g, String.format("M%02d_Lagrangian_u_g%02d.dat",s_i,this.iterNum));
		MathFunc fa = new Vector2MathFunc(a);

		
		if(bTestWholdDomain)
			weakForm.setF(fu_g.M(-1.0));//!!!!!!1.0=>-1.0???   //???.A(this.modelReal.delta)
		else
			weakForm.setF(FC.c(0.0));

		
		//d*u + k*u_n = q
		//采用自然边界：u_n + u = 0
		if(bTestWholdDomain) {
			weakForm.setParam(
					FC.C1.D(fa.M(model_k)),
					FC.C1,
					null,//q=0
					//FC.c0
					//***
					FC.C1.D(fa.M(model_k))//d=k
				);
		} else {
			//d*u + k*u_n = q
			//自然边界(u_n+u=0)+边界测量(-(u-g,\psi)_\Gamma)
			//
			weakForm.setParam(
					FC.C1.D(fa.M(model_k)),
					FC.C1,
					fu_g.M(-1.0), //q=-(u-g)
					FC.C1.D(fa.M(model_k)) //d=k
				);
		}
		mesh.clearBorderNodeMark();
		HashMap<NodeType, MathFunc> mapNTF = new HashMap<NodeType, MathFunc>();
		if(bTestWholdDomain && bTestWholeDomainDirichletBoundary)
			mapNTF.put(NodeType.Dirichlet, null);
		else
			mapNTF.put(NodeType.Robin, null);
		mesh.markBorderNode(mapNTF);
		
		AssemblerScalar assembler = new AssemblerScalar(mesh, weakForm);
		System.out.println("Begin Assemble...lambda");
		assembler.assemble();
		SparseMatrix stiff = assembler.getStiffnessMatrix();
		SparseVector load = assembler.getLoadVector();
		if(bTestWholdDomain && bTestWholeDomainDirichletBoundary)
			assembler.imposeDirichletCondition(FC.C0);
		System.out.println("Assemble done!");

		Equation eqn = new Equation();
		eqn.A = stiff;
		eqn.f = load;
		
		return eqn;
	}	
	
	/**
	 * Get stiff matrix and load vector for equation of 'u'
	 * 
	 * L_{\lambda}(u)=0:
	 *   ((1/(a*k))*\nabla{u},\nabla{\phi}) + (u,\phi) = 0
	 *  
	 * Boundary Condition: Robin
	 *   (1/(a*k))*\partial_{n}{u}  + (1/(a*k))*u = 0, on \Gamma
	 *   =>（实际计算的时候不能，参数中不能消去1/(a*k)）
	 *   (1/(a*k))*(\partial_{n}{u}  + u) = 0, on \Gamma
	 * 
	 * @param a
	 * @param g: 可选项
	 *     g==null, 以自然边界条件在大区域上计算u，然后截取到小区域上
	 *     g!=null, 以g为Dirichlet边界条件在小区域上计算u
	 * @return
	 */
	public void getOrSolveEqnU(
			Mesh _mesh, Vector a, Vector g, Vector u0_x, Vector u0_y, //In
			Equation eqn, Vector u) { //Out
			 
		
		WeakFormLaplace2D weakForm = new WeakFormLaplace2D();
		MathFunc fa = new Vector2MathFunc(a);
		
		//不能忽略光源的影响???this.bTestWholeDomainDirichletBoundary
		if(this.useVectorMu_a) {//2011/10/18
			Vector aRealNew = Tools.interplateFrom(aMesh, _mesh, 
					new Vector2MathFunc(aReal,aMesh,"x","y").
					setDefaultFunction(FC.c(this.aBackground)));
			modelReal.setMu_a(new Vector2MathFunc(aRealNew));
		}
		//if(g == null)
			weakForm.setF(this.modelReal.getDelta());
		//else
		//	weakForm.setF(FC.c(0.0));
		
		//DuDn du0dn = new DuDn(new Vector2Function(u0_x),new Vector2Function(u0_y),null);
		
		weakForm.setParam(
				FC.C1.D(fa.M(model_k)),
				FC.C1,
				null,
				FC.C1.D(fa.M(model_k)));
		
//		weakForm.setParam(
//				FC.c1.D(fa.M(model_k)),
//				FC.c1,
//				FC.c1.D(fa.M(model_k)).M(du0dn),
//				FC.c0);
		
		_mesh.clearBorderNodeMark();
		HashMap<NodeType, MathFunc> mapNTF = new HashMap<NodeType, MathFunc>();
		if(g == null)
			mapNTF.put(NodeType.Robin, null);
		else
			mapNTF.put(NodeType.Dirichlet, null);
		_mesh.markBorderNode(mapNTF);

		AssemblerScalar assembler = new AssemblerScalar(_mesh, weakForm);
		System.out.println("Begin Assemble...u");
		assembler.assemble();
		SparseMatrix stiff = assembler.getStiffnessMatrix();
		SparseVector load = assembler.getLoadVector();
		if(g != null)
			assembler.imposeDirichletCondition(new Vector2MathFunc(g));
		System.out.println("Assemble done!");

		if(eqn != null) {
			eqn.A = stiff;
			eqn.f = load;
		}
		if(u != null) {
	        SolverJBLAS sol = new SolverJBLAS();
			Vector uSol = sol.solveDGESV(stiff, load);
			u.set(uSol);
		}
	}
	
	public Equation getEqnU(Vector a, Vector g, Vector u0_x, Vector u0_y) {
		if(g == null) {
			//Robin条件，由于光源在区域外面，在小区域上直接求解会得到0解，因此
			//先在大区域上求解uBig，然后截取解到小区域uSmall，
			//最后将uSmall在边界上的值作为Dirichlet边界在小区域上求解
			Vector aBig = Tools.extendData(mesh, meshBig, a, this.aBackground);
			if(debug)
				plotVector(meshBig,aBig,String.format("aBig_ext%02d.dat",this.iterNum));
			
	        Vector uBig = new SparseVectorHashMap();
			getOrSolveEqnU(meshBig, aBig, null, null, null, null, uBig);
			
			if(debug)
				plotVector(meshBig,uBig,String.format("uBig_ext%02d.dat",this.iterNum));
	        Vector uSmall = Tools.extractData(meshBig, mesh, uBig);
	        
	        Equation eq = new Equation();
	        getOrSolveEqnU(mesh,a,uSmall, null, null, eq, null);
	        //getEqnU(mesh,a,null, u0_x, u0_y, eq, null);
	        return eq;
		}
		else {
			//在小区域求解Dirichlet问题（得到系数矩阵和右端向量）
	        Equation eq = new Equation();
			getOrSolveEqnU(mesh,a,g, null, null, eq, null);
			return eq;
		}
	}
	
	
	/**
	 * Get stiff matrix and load vector for equation of 'a(x)'
	 * 
	 * L_a(a)=0:
	 * 
	 * \beta(a-aGlob,\chi) - 
	 * 					\sum_{i=1,N}{
	 * 						 ((ak*k)^{-2}*k*\chi\nabla{u_i},\nabla{\lambda_i})
	 * 					} = 0
	 * =>
	 * (a,\chi) = (aGlob,\chi) + (1/\beta)*
	 *                  \sum_{i=1,N}{
	 * 						 ((ak*k)^{-2}*k*\chi*\nabla{u_i},\nabla{\lambda_i})
	 * 					}
	 */
	public Equation getEqnA(Vector[] u, Vector[] u_x,Vector[] u_y, Vector[] lambda, 
			Vector ak, Vector _aGlob, 
			MathFunc diri) {
//        WeakFormLaplace2D weakForm = new WeakFormLaplace2D();
//        weakForm.setParam(FC.c0, FC.c1, null, null);
//        //stabilize
//        //weakForm.setParam(FC.c(0.001), FC.c1, null, null);
//        
//        Function faGlob = new Vector2Function(_aGlob);
//        int N=u.length;
//        Vector[] uDotlmd = new Vector[N];
//        for(int i=0; i<N; i++) {
//        	Vector ux = Tools.computeDerivative(mesh, u[i], "x");
//        	Vector uy = Tools.computeDerivative(mesh, u[i], "y");
//        	Vector lx = Tools.computeDerivative(mesh, lambda[i], "x");
//        	Vector ly = Tools.computeDerivative(mesh, lambda[i], "y");
//        	plotVector(mesh, ux, String.format("M%02d_La_RHS_ux%02d.dat",i,this.iterNum));
//        	plotVector(mesh, uy, String.format("M%02d_La_RHS_uy%02d.dat",i,this.iterNum));
//        	plotVector(mesh, lx, String.format("M%02d_La_RHS_lx%02d.dat",i,this.iterNum));
//        	plotVector(mesh, ly, String.format("M%02d_La_RHS_ly%02d.dat",i,this.iterNum));
//        	uDotlmd[i] = new SparseVector(ux.getDim());
//        	uDotlmd[i].add(FMath.axMuly(1.0, ux, lx));
//        	uDotlmd[i].add(FMath.axMuly(1.0, uy, ly));
//        	plotVector(mesh,uDotlmd[i],
//        			String.format("M%02d_La_RHS_v_lmd_Grad%02d.dat",i,this.iterNum));
//        	//this.connectSells(mesh, uDotlmd[i]);
//        	//plotVector(mesh,uDotlmd[i],
//        	//		String.format("M%02d_La_RHS_v_lmd_Grad_ConnectSells%02d.dat",i,this.iterNum));
//        }
//		Vector sum2 = FMath.sum(uDotlmd);
//		//this.connectSells(mesh, sum2);
//		plotVector(mesh,sum2,String.format("La_RHS_v_lmd_sum2_%02d.dat",this.iterNum));
//
//		Function akmk_2mk = FMath.pow(new Vector2Function(ak).M(model_k),-2.0).M(model_k);
//		plotFunction(mesh,akmk_2mk,String.format("La_RHS_akmk_2_%02d.dat",this.iterNum));
//		
//		Function rhs = akmk_2mk.M(new Vector2Function(sum2)).M(1.0/beta);
//		plotFunction(mesh,rhs,String.format("La_RHS_rhs%02d.dat",this.iterNum));
//		Function f2 = faGlob.A(rhs);
//		plotFunction(mesh,f2,String.format("La_RHS_all%02d.dat",this.iterNum));
//		weakForm.setF(f2);
		
		MathFunc faGlob = new Vector2MathFunc(_aGlob);
		MathFunc akmk_2mk = FMath.pow(new Vector2MathFunc(ak).M(model_k),-2.0)
							.M(model_k).M(1.0/beta);
		//Function ak_2 = FMath.pow(new Vector2Function(ak),-2.0).M(1.0/beta);
		int NF = u.length;
		MathFunc[] fu = new MathFunc[NF];
		MathFunc[] fl = new MathFunc[NF];
		MathFunc[] fumfl = new MathFunc[NF];
		DuDn[] du0dn = new DuDn[NF];
		MathFunc[] du0dnmfl = new MathFunc[NF];
		for(int k=0;k<NF;k++) {
			fu[k] = new Vector2MathFunc(u[k]);
			fl[k] = new Vector2MathFunc(lambda[k]);
			fumfl[k] = fu[k].M(fl[k]);
        	du0dn[k] = new DuDn(
        		new Vector2MathFunc(u_x[k]),
        		new Vector2MathFunc(u_y[k]),
        		null
        		);
        	du0dnmfl[k] = du0dn[k].M(fl[k]);
        }
		
        WeakFormLa weakForm = new WeakFormLa();
        
        //***
        //weakForm.setRobin(akmk_2mk.M(FMath.sum(fumfl)), FC.c0);
        //plotFunction(mesh, akmk_2mk.M(FMath.sum(fumfl)), "La_robin.dat");
        weakForm.setRobin(FC.C0, FC.C0);
        //weakForm.setRobin(akmk_2mk.M(FMath.sum(du0dnmfl)).M(-1.0), FC.c0);
        
        weakForm.setF(faGlob,//FC.c0,//.A(ak_2.M(FMath.sum(fl)).M(modelReal.delta)) 
        		akmk_2mk, FC.C0,
        		fu, fl
        	);
		

        //需要重新标记边界条件
		mesh.clearBorderNodeMark();
		HashMap<NodeType, MathFunc> mapNTF = new HashMap<NodeType, MathFunc>();
		if(this.bTestWholeDomainDirichletBoundary)
			mapNTF.put(NodeType.Dirichlet, null);
		else
			mapNTF.put(NodeType.Robin, null);
		mesh.markBorderNode(mapNTF);
		
        //5.Assembly process
        AssemblerScalar assembler =
                new AssemblerScalar(mesh, weakForm);
        System.out.println("Begin Assemble...a(x)");
        assembler.assemble();
        SparseMatrix stiff = assembler.getStiffnessMatrix();
        SparseVector load = assembler.getLoadVector();
        
        //Boundary condition
        if(this.bTestWholeDomainDirichletBoundary)
        	assembler.imposeDirichletCondition(diri);
        
        System.out.println("Assemble done!");

		Equation eqn = new Equation();
		eqn.A = stiff;
		eqn.f = load;
		
		return eqn;
	}
	
//	public Equation getEqnA(Vector u, Vector lambda, 
//			Vector ak, Vector _aGlob, 
//			Function diri) {
//		Vector[] vu = new Vector[1];
//		Vector[] vlmd = new Vector[1];
//		vu[0] = u;
//		vlmd[0] = lambda;
//		return getEqnA(vu, vlmd, ak, _aGlob, diri);
//	}
	
	
	/**
	 * Residual of state equation L_{\lambda}(u)
	 * 
	 * 获取状态方程的余量
	 * @return
	 */
	public SparseVector getResLlmd(Vector a, Vector u, Vector g, Vector u0_x, Vector u0_y) {
		Equation eq = this.getEqnU(a, g, u0_x, u0_y);
        SparseVector res = eq.f.copy();
        res.setAll(0.0);
        
        eq.A.mult(u, res);
        res.add(-1.0, eq.f);
        //connectSells(mesh, res);
        zeroHangingNode(mesh,res);
        return res;
	}
	
	/**
	 * Residual of adjoint equation L_u(\lambda)
	 * 
	 * @return
	 */
	public SparseVector getResLu(int s_i,
			Vector a, Vector u, Vector g,
			Vector lambda) {
		Equation eq = this.getEqnLambda(s_i, a, u, g);
        SparseVector res = eq.f.copy();
        res.setAll(0.0);
        
        eq.A.mult(lambda, res);
        res.add(-1.0, eq.f);
        //res.axpy(-1.0, eq.f);
        //connectSells(mesh, res);
        zeroHangingNode(mesh,res);
        return res;
    }
	
	/**
	 * Residual of parameter regularization L_{q}
	 * 
	 * @param v
	 * @param lambda
	 * @param a_glob
	 * @return
	 */
//	public Vector getResLq(Vector u, Vector lambda,
//			Vector ak, Vector _aGlob, Function diri) {
//		Equation eq = this.getEqnA(u, lambda, ak, _aGlob, diri);
//        Vector res = new SparseVector(eq.f.getDim());
//        eq.A.mult(ak, res);
//        res.add(-1.0, eq.f);
//        //connectSells(mesh, res);
//        //zeroHangingNode(mesh,res);
//        plotVector(mesh,res,String.format("Res_La%02d.dat",this.iterNum));
//        
//        //光滑余量
//        //res = Utils.gaussSmooth(mesh, res, 2, 0.5);
//        //plotVector(mesh,res,String.format("Res_LaSmooth%02d.dat",this.iterNum));
//        
//        //Solver sol = new Solver();
//        //Vector a_solve = sol.solveCGS(eq.A, eq.f);
//        SolverJBLAS sol = new SolverJBLAS();
//		Vector a_solve = sol.solveDGESV(eq.A, eq.f);
//		
//        plotVector(mesh,a_solve,String.format("a_solve%02d.dat",this.iterNum));
//        return res;
//	}
	
	public SparseVector getResLq(Vector[] u,Vector[] u_x,Vector[] u_y, Vector[] lambda,
			Vector ak, Vector _aGlob, MathFunc diri) {
		Equation eq = this.getEqnA(u, u_x, u_y, lambda, ak, _aGlob, diri);
        SparseVector res = eq.f.copy();
        res.setAll(0.0);
        
        eq.A.mult(ak, res);
        plotVector(mesh,res,String.format("Res_La_mult%02d.dat",this.iterNum));
        //this.connectSells(mesh, res);
        //plotVector(mesh,res,String.format("Res_La_mult_connectSells%02d.dat",this.iterNum));
        
        plotVector(mesh,eq.f,String.format("Res_La_f%02d.dat",this.iterNum));
        //this.connectSells(mesh, eq.f);
        //plotVector(mesh,eq.f,String.format("Res_La_f_connectSells%02d.dat",this.iterNum));

        res.add(-1.0, eq.f);
        plotVector(mesh,res,String.format("Res_La%02d.dat",this.iterNum));
        //this.connectSells(mesh, res);
        //plotVector(mesh,res,String.format("Res_La_connectSells%02d.dat",this.iterNum));
        
        zeroHangingNode(mesh,res);
        plotVector(mesh,res,String.format("Res_La_zeroHangingNode%02d.dat",this.iterNum));

        
        //光滑余量
        //res = Utils.gaussSmooth(mesh, res, 1, 0.5);
        //plotVector(mesh,res,String.format("Res_LaSmooth%02d.dat",this.iterNum));
        
        //直接求解a(x)
        //Solver sol = new Solver();
        //Vector a_solve = sol.solveCGS(eq.A, eq.f);
        SolverJBLAS sol = new SolverJBLAS();
		Vector a_solve = sol.solveDGESV(eq.A, eq.f);
		
        plotVector(mesh,a_solve,String.format("a_solve%02d.dat",this.iterNum));
        return res;
	}
	
	
	/////////////////////////////////////////////////////////////////////////
	//Begin*************Construction of search direction matrix**************
	/*
	 *  ( M  A'  0 )(du)     (Lu)
	 *  ( A  0   C )(dl) = - (Ll)
	 *  ( 0  C' bR )(dq)     (Lq)
	 */
	//***********************************************************************
	
	/**
	 * Weak form of 'M'
	 * 
	 * 整体测量：(du,\phi)_\Omega
	 * 边界测量：(du,\phi)_\Gamma
	 * 
	 */
	public SparseMatrix getM() {
        WeakFormLaplace2D weakForm = new WeakFormLaplace2D();
        
        if(this.bTestWholdDomain) {
            //weakForm.setParam(FC.c0, FC.c1, null, null);
            //stabilize 使结果更光滑一些，不会导致结果有大的变化
            weakForm.setParam(
            		FC.c(0.0), //---同时改
            		FC.c(1.0), 
            		null,
            		FC.c(0.0)); //---同时改
        } else {
        	weakForm.setParam(
        		FC.c(0.0), //null==(k=1)
        		FC.c(0.0), 
        		null, //q=0
        		FC.C1 //d=1
        		);
        }

        //Right hand side(RHS): f(x) = 0
        weakForm.setF(FC.C0);

        //需要重新标记边界条件，否则在“整体合成”过程会出现错误。
        //虽然边界条件实在大矩阵中设置，这一步也是需要的。
		mesh.clearBorderNodeMark();
		HashMap<NodeType, MathFunc> mapNTF = new HashMap<NodeType, MathFunc>();
		
        if(this.bTestWholdDomain) 
        	mapNTF.put(NodeType.Dirichlet, null);
		else
			mapNTF.put(NodeType.Robin, null);
		
		mesh.markBorderNode(mapNTF);

        //5.Assembly process
        AssemblerScalar assembler =
                new AssemblerScalar(mesh, weakForm);
        System.out.println("Begin Assemble...M");
        assembler.assemble(false);
        SparseMatrix stiff = assembler.getStiffnessMatrix();
        //Boundary condition 
        //边界条件需要在大矩阵中设置
        //assembler.imposeDirichletCondition(FC.c0);
        System.out.println("Assemble done!");
        
        //正则化 \eps*I + M
//        System.out.println("M");
//        for(int i=1;i<=stiff.getColDim();i++) {
//        	System.out.println(stiff.get(i, i));
//        	stiff.set(i, i, stiff.get(i, i)*1.9);
//        }
        return stiff;
    }
	
	/**
	 * Weak form of 'A'
	 * 
	 * ((1/(a*k)*\nabla{du},\nabla{\psi}) + (du,\phi)
	 * 
	 * where 
	 *   du=\delta{u}
	 */
	public SparseMatrix getA(Vector ak, boolean procHangingNode) {
		return getA(ak,null,null,procHangingNode).A;
	}
	public Equation getA(Vector ak, MathFunc f, MathFunc diri,boolean procHangingNode) {
		WeakFormLaplace2D weakForm = new WeakFormLaplace2D();
		MathFunc fa = new Vector2MathFunc(ak);
		
		if(f==null)
			weakForm.setF(FC.c(0.0));
			//weakForm.setF(this.modelReal.delta);
		else
			weakForm.setF(f);
		
		weakForm.setParam(
				FC.C1.D(fa.M(model_k)),
				FC.C1,
				null,
				//***
				FC.C1.D(fa.M(model_k))
				//FC.c0
				);


        //需要重新标记边界条件，否则在“整体合成”过程会出现错误。
        //虽然边界条件实在大矩阵中设置，这一步也是需要的。
		mesh.clearBorderNodeMark();
		HashMap<NodeType, MathFunc> mapNTF = new HashMap<NodeType, MathFunc>();
		
		if(this.bTestWholeDomainDirichletBoundary)
			mapNTF.put(NodeType.Dirichlet, null);
		else
			mapNTF.put(NodeType.Robin, null);//bugfix add
		mesh.markBorderNode(mapNTF);
		
        //5.Assembly process
        AssemblerScalar assembler =
                new AssemblerScalar(mesh, weakForm);
        System.out.println("Begin Assemble...A");
        assembler.assemble(procHangingNode);
        SparseMatrix stiff = assembler.getStiffnessMatrix();
        SparseVector load = null;
        if(f != null)
        	load = assembler.getLoadVector();
        if(diri != null)//'A'的边界条件需要在大矩阵中设置
        	assembler.imposeDirichletCondition(diri);
        System.out.println("Assemble done!");
        
		Equation eqn = new Equation();
		eqn.A = stiff;
		eqn.f = load;
		
		//stiff.print();
//        //
//		System.out.println("A");
//        for(int i=1;i<=stiff.getColDim();i++) {
//        	System.out.println(stiff.get(i, i));
//        }
        return eqn;
	}
	
	/**
	 * Weak form of 'AT'
	 * 
	 * ((1/(a+k)*\nabla{\phi},\nabla{dl}) + ((a*\phi,dl)
	 * 
	 * where 
	 *   dl=\delta{\lambda}
	 */
	public SparseMatrix getAT(Vector ak) {
		return getA(ak,true); //Note for adaptive mesh hanging nodes
		//return getA(ak,true).trans();
	}
	
	/**
	 * Weak form of 'C'
	 * 
	 * ((-(a*k)^{-2}*k*da*\nabla{u},\nabla{\psi})
	 * 
	 * @param ak
	 * @param uk
	 * @return
	 */
	public SparseMatrix getC(Vector ak, Vector uk, Vector uk_x, Vector uk_y) {
//        WeakFormGCMDual weakForm = new WeakFormGCMDual();
//        Function fa = new Vector2Function(ak);
//        Function b1 = new Vector2Function(Tools.computeDerivative(mesh, uk, "x"));
//        Function b2 = new Vector2Function(Tools.computeDerivative(mesh, uk, "y"));
//        Function _amk_2mk = FMath.pow(fa.M(model_k),-2).M(model_k).M(-1.0);
//        
//        weakForm.setParam(FC.c0, FC.c0, 
//        		b1.M(_amk_2mk), 
//        		b2.M(_amk_2mk));
        
        WeakFormC weakForm = new WeakFormC();
        MathFunc fa = new Vector2MathFunc(ak);
        MathFunc fu = new Vector2MathFunc(uk);
        MathFunc _amk_2mk = FMath.pow(fa.M(model_k),-2).M(model_k).M(-1.0);
        
        weakForm.setParam(_amk_2mk, FC.C0, fu);
        //weakForm.setParam(_amk_2mk, fu, fu);
        
        //***
        //Robin:  d*u + k*u_n= q (自然边界：d==k, q=0)
        //weakForm.setRobin(FC.c0, _amk_2mk.M(fu));
        //DuDn du0dn = new DuDn(
        //		new Vector2Function(uk_x),
        //		new Vector2Function(uk_y),
        //		null);
        //weakForm.setRobin(FC.c0, _amk_2mk.M(du0dn).M(-1.0));
        weakForm.setRobin(FC.C0, FC.C0);
        
        //stabilize
        //weakForm.setParam(FC.c(0.0001), FC.c0, b1.M(_apk_2), b2.M(_apk_2));
        weakForm.setF(FC.C0);

        //需要重新标记边界条件，否则在“整体合成”过程会出现错误。
        //虽然边界条件实在大矩阵中设置，这一步也是需要的。
		mesh.clearBorderNodeMark();
		HashMap<NodeType, MathFunc> mapNTF = new HashMap<NodeType, MathFunc>();
		if(this.bTestWholeDomainDirichletBoundary)
			mapNTF.put(NodeType.Dirichlet, null);
		else
			mapNTF.put(NodeType.Robin, null);
		mesh.markBorderNode(mapNTF);
		
        //5.Assembly process
        AssemblerScalar assembler =
                new AssemblerScalar(mesh, weakForm);
        System.out.println("Begin Assemble...C");
        assembler.assemble(false);
        SparseMatrix stiff = assembler.getStiffnessMatrix();
        //Boundary condition
        //边界条件需要在大矩阵中设置
        //assembler.imposeDirichletCondition(FC.c0);
        System.out.println("Assemble done!");
        
        //stiff.print();
        
//        System.out.println("C");
//        for(int i=1;i<=stiff.getColDim();i++) {
//        	System.out.println(stiff.get(i, i));
//        }
        return stiff;
	}
	
	/**
	 * Weak form of 'CT'
	 * 
	 * ((-(a*k)^{-2}*k*\chi*\nabla{u},\nabla{dl})
	 * 
	 * where
	 * 
	 *   dl=\delta{\lambda}
	 * 
	 * @param ak
	 * @param uk
	 * @return
	 */
	public SparseMatrix getCT(Vector ak, Vector uk, Vector uk_x, Vector uk_y) {
//        WeakFormGCM weakForm = new WeakFormGCM();
//        Function fa = new Vector2Function(ak);
//        Function b1 = new Vector2Function(Tools.computeDerivative(mesh, uk, "x"));
//        Function b2 = new Vector2Function(Tools.computeDerivative(mesh, uk, "y"));
//        Function _amk_2mk = FMath.pow(fa.M(model_k),-2).M(model_k).M(-1.0);
//        
//        weakForm.setParam(FC.c0, FC.c0, 
//        		b1.M(_amk_2mk), 
//        		b2.M(_amk_2mk));
        
        WeakFormCT weakForm = new WeakFormCT();
        MathFunc fa = new Vector2MathFunc(ak);
        MathFunc fu = new Vector2MathFunc(uk);
        MathFunc _amk_2mk = FMath.pow(fa.M(model_k),-2).M(model_k).M(-1.0);
        
        weakForm.setParam(_amk_2mk, FC.C0, fu);
        //weakForm.setParam(_amk_2mk, fu, fu);
        
        //***
        //Robin:  d*u + k*u_n= q (自然边界：d==k, q=0)
        //weakForm.setRobin(FC.c0, _amk_2mk.M(fu));
        //DuDn du0dn = new DuDn(
        //		new Vector2Function(uk_x),
        //		new Vector2Function(uk_y),
        //		null);
        //2011/1/18 这两个条件结果差不多（bugfix:是因为没有标记边界条件：忘记调用mesh.markBorderNode(mapNTF);）
        //weakForm.setRobin(FC.c0, _amk_2mk.M(du0dn).M(-1.0));
        weakForm.setRobin(FC.C0, FC.C0);
        
        //stabilize
        //weakForm.setParam(FC.c(0.0001), FC.c0, b1.M(_apk_2), b2.M(_apk_2));
        weakForm.setF(FC.C0);

        //需要重新标记边界条件，否则在“整体合成”过程会出现错误。
        //虽然边界条件实在大矩阵中设置，这一步也是需要的。
		mesh.clearBorderNodeMark();
		HashMap<NodeType, MathFunc> mapNTF = new HashMap<NodeType, MathFunc>();
		if(this.bTestWholeDomainDirichletBoundary)
			mapNTF.put(NodeType.Dirichlet, null);
		else
			mapNTF.put(NodeType.Robin, null);
		mesh.markBorderNode(mapNTF);
		
        //5.Assembly process
        AssemblerScalar assembler =
                new AssemblerScalar(mesh, weakForm);
        System.out.println("Begin Assemble...CT");
        assembler.assemble(false);
        SparseMatrix stiff = assembler.getStiffnessMatrix();
        //Boundary condition
        //边界条件需要在大矩阵中设置
        //assembler.imposeDirichletCondition(FC.c0);
        System.out.println("Assemble done!");
        
        //stiff.print();
        
        return stiff;		
	}
	
	public Matrix testGetCT(Vector ak, Vector uk, Vector uk_x, Vector uk_y) {
      WeakFormCT weakForm = new WeakFormCT();
      MathFunc fa = new Vector2MathFunc(ak);
      MathFunc fu = new Vector2MathFunc(uk);
      MathFunc _amk_2mk = FMath.pow(fa.M(model_k),-2).M(model_k).M(-1.0);
      
      weakForm.setParam(_amk_2mk, fu, fu);
      //***
      //Robin:  d*u + k*u_n= q (自然边界：d==k, q=0)
      //weakForm.setRobin(FC.c0, _amk_2mk.M(fu));
      //DuDn du0dn = new DuDn(
      //		new Vector2Function(uk_x),
      //		new Vector2Function(uk_y),
      //		null);
      //2011/1/18 这两个条件结果差不多（bugfix:是因为没有标记边界条件：忘记调用mesh.markBorderNode(mapNTF);）
      //weakForm.setRobin(FC.c0, _amk_2mk.M(du0dn).M(-1.0));
      weakForm.setRobin(FC.C0, FC.C0);

      //stabilize
      //weakForm.setParam(FC.c(0.0001), FC.c0, b1.M(_apk_2), b2.M(_apk_2));
	  weakForm.setF(FC.C0);

      
      //需要重新标记边界条件，否则在“整体合成”过程会出现错误。
      //虽然边界条件实在大矩阵中设置，这一步也是需要的。
		mesh.clearBorderNodeMark();
		HashMap<NodeType, MathFunc> mapNTF = new HashMap<NodeType, MathFunc>();
		if(this.bTestWholeDomainDirichletBoundary)
			mapNTF.put(NodeType.Dirichlet, null);
		else
			mapNTF.put(NodeType.Robin, null);
		mesh.markBorderNode(mapNTF);
		
      //5.Assembly process
      AssemblerScalar assembler =
              new AssemblerScalar(mesh, weakForm);
      System.out.println("Begin Assemble...CT");
      assembler.assemble(false);
      Matrix stiff = assembler.getStiffnessMatrix();
      //Boundary condition
      //边界条件需要在大矩阵中设置
      //assembler.imposeDirichletCondition(FC.c0);
      System.out.println("Assemble done!");
      
      //stiff.print();
      
      return stiff;		
	}	
	public Equation testGetCTLoad(Vector[] u, Vector[] u_x,Vector[] u_y, Vector[] lambda, 
			Vector ak, MathFunc diri) {
		
		MathFunc akmk_2mk = FMath.pow(new Vector2MathFunc(ak).M(model_k),-2.0).M(model_k).M(-1.0/beta);
		int NF = u.length;
		MathFunc[] fu = new MathFunc[NF];
		MathFunc[] fl = new MathFunc[NF];
		for(int k=0;k<NF;k++) {
			fu[k] = new Vector2MathFunc(u[k]);
			fl[k] = new Vector2MathFunc(lambda[k]);
        }
        WeakFormLa weakForm = new WeakFormLa();
        
        weakForm.setRobin(FC.C0, FC.C0);
        weakForm.setF(FC.C0,
        		akmk_2mk, FC.C1,
        		fu, fl
        	);
		

        //需要重新标记边界条件
		mesh.clearBorderNodeMark();
		HashMap<NodeType, MathFunc> mapNTF = new HashMap<NodeType, MathFunc>();
		if(this.bTestWholeDomainDirichletBoundary)
			mapNTF.put(NodeType.Dirichlet, null);
		else
			mapNTF.put(NodeType.Robin, null);
		mesh.markBorderNode(mapNTF);
		
        //5.Assembly process
        AssemblerScalar assembler =
                new AssemblerScalar(mesh, weakForm);
        System.out.println("Begin Assemble...a(x)");
        assembler.assemble();
        SparseMatrix stiff = assembler.getStiffnessMatrix();
        SparseVector load = assembler.getLoadVector();
        
        //Boundary condition
        if(this.bTestWholeDomainDirichletBoundary)
        	assembler.imposeDirichletCondition(diri);
        
        System.out.println("Assemble done!");

		Equation eqn = new Equation();
		eqn.A = stiff;
		eqn.f = load;
		
		return eqn;
	}
	
	
	/**
	 * Weak form of '\beta*R'
	 * 
	 * \beta*(da,\chi)

	 * @return
	 */
	public SparseMatrix getBR() {
        WeakFormLaplace2D weakForm = new WeakFormLaplace2D();
        weakForm.setParam(FC.C0, FC.c(beta), null, null);
        //stabilize
        //weakForm.setParam(FC.c(1000), FC.c(beta), null, null);
        
        weakForm.setF(FC.C0);

        //需要重新标记边界条件，否则在“整体合成”过程会出现错误。
        //虽然边界条件实在大矩阵中设置，这一步也是需要的。
		mesh.clearBorderNodeMark();
		HashMap<NodeType, MathFunc> mapNTF = new HashMap<NodeType, MathFunc>();
		mapNTF.put(NodeType.Dirichlet, null);
		mesh.markBorderNode(mapNTF);
		
        //5.Assembly process
        AssemblerScalar assembler =
                new AssemblerScalar(mesh, weakForm);
        System.out.println("Begin Assemble...R");
        assembler.assemble();
        SparseMatrix stiff = assembler.getStiffnessMatrix();
        //Boundary condition
        //边界条件需要在大矩阵中设置
        //assembler.imposeDirichletCondition(FC.c0);
        System.out.println("Assemble done!");
        
//        //
//        System.out.println("\beta R");
//        for(int i=1;i<=stiff.getColDim();i++) {
//        	System.out.println(stiff.get(i, i));
//        }
        return stiff;
	}
	
	//End*************Construction of search direction matrix**************
	/////////////////////////////////////////////////////////////////////////
	
	/**
	 * 真解u
	 * 
	 * @return
	 */
	public Vector solveRealU(int s_i) {
		if(this.useVectorMu_a) {//2011/10/18
			Vector aRealBig = Tools.interplateFrom(aMesh, meshBig, 
					new Vector2MathFunc(aReal,aMesh,"x","y").
					setDefaultFunction(FC.c(this.aBackground)));
			modelReal.setMu_a(new Vector2MathFunc(aRealBig));
		}
		Vector uRealBig = modelReal.solveNeumann(meshBig);
		plotVector(meshBig, uRealBig, String.format("M%02d_uRealBig.dat",s_i));
		
		//截取meshBig的部分解到mesh上
		Vector uReal = Tools.extractData(meshBig, mesh, uRealBig);
		plotVector(mesh, uReal, String.format("M%02d_uReal.dat",s_i));

//以下验证都成功（2011/8/4）		
//		//验证从大区域截取出来的解与边界施加Dirichlet条件解是否相同
//		Vector aRealVec = Tools.function2vector(mesh, modelReal.mu_a);
//		Vector uRealDiri = solveStateEquation(aRealVec, uReal);
//		plotVector(mesh, uRealDiri, String.format("M%02d_uRealDiri.dat",s_i));
//		plotVector(mesh, FMath.axpy(-1.0, uRealDiri,uReal), String.format("M%02d_uReal_uRealDiri_diff.dat",s_i));
//		
//		//比较解：边界相同uSmall，a(x)不同
//		Vector aGuessVec = Tools.function2vector(mesh, modelGuess.mu_a);
//		Vector uGuessDiriReal = solveStateEquation(aGuessVec, uReal);
//		Vector uGuessDiriReal2 = modelGuess.solveDirichlet(mesh, new Vector2Function(uReal));
//		plotVector(mesh, uGuessDiriReal, String.format("M%02d_uGuessDiriReal.dat",s_i));
//		plotVector(mesh, FMath.axpy(-1.0, uGuessDiriReal,uRealDiri), String.format("M%02d_uGuessDiriReal_uRealDiri_diff.dat",s_i));
//		plotVector(mesh, FMath.axpy(-1.0, uGuessDiriReal2,uRealDiri), String.format("M%02d_uGuessDiriReal2_uRealDiri_diff.dat",s_i));
//		
//		//比较解：a(x)相同aGuessVec，边界不同
//		Vector uGuessBig = modelGuess.solveNeumann(meshBig);
//		plotVector(meshBig, uGuessBig, String.format("M%02d_uGuessBig.dat",s_i));
//		Vector uGuess = Tools.extractData(meshBig, mesh, uGuessBig);
//		plotVector(mesh, uGuess, String.format("M%02d_uGuess.dat",s_i));
//		plotVector(mesh, FMath.axpy(-1.0, uGuessDiriReal,uGuess), String.format("M%02d_uGuessDiriReal_uGuess_diff.dat",s_i));
		
        return uReal;
	}
	
	/**
	 * Solve du=\delat{u} based on the second search direction equation
	 * 
	 *     A*du + C*da = -ResLlmd
	 *   =>
	 *     A*du = -ResLlmd-C*da
	 *   =>
	 *     du = inv(A)*(-ResLlmd-C*da)
	 *     
	 * where 
	 *  ResLlmd = residual of L_{\lambda}
	 * 
	 * @param ak
	 * @param _resLlmd_da: -ResLlmd-C*\delta{a}
	 * @param uk
	 * @return
	 */
	public Vector solveDeltaU(Vector ak, Vector _resLlmd_da, Vector uk) {
		Equation eq = getA(ak,new Vector2MathFunc(_resLlmd_da),FC.C0,true);
        //Solver sol = new Solver();
        //Vector x = sol.solveCGS(eq.A, eq.f);
        SolverJBLAS sol = new SolverJBLAS();
		Vector x = sol.solveDGESV(eq.A, eq.f);
        return x;
	}

	/**
	 * 求解关于u的状态方程
	 * 
	 * @param a 
	 * @return
	 */
	public Vector solveStateEquation(Vector a, Vector g, Vector u0_x, Vector u0_y) {
		Vector u = new SparseVectorHashMap();
		if(g == null) {
			//Robin条件，由于光源在区域外面，在小区域上直接求解会得到0解，因此
			//先在大区域上求解，然后截取解到小区域
			Vector aBig = Tools.extendData(mesh, meshBig, a, this.aBackground);
			if(debug)
				plotVector(meshBig,aBig,String.format("aBig_ext%02d.dat",this.iterNum));
			
	        Vector uBig = new SparseVectorHashMap();
			getOrSolveEqnU(meshBig,aBig,null, null, null, null, uBig);
			if(debug)
				plotVector(meshBig,uBig,String.format("uBig_ext%02d.dat",this.iterNum));
	        u = Tools.extractData(meshBig, mesh, uBig);
		} else {
			getOrSolveEqnU(mesh, a, g, u0_x, u0_y, null, u);
		}
		return u;
	}
	
	/**
	 * 求解关于lambda的伴随方程
	 * 
	 * @param s_i
	 * @param a
	 * @param u
	 * @param g: u|_\Gamma
	 * @return
	 */
	public Vector solveAdjointEquation(int s_i,  
			Vector a, Vector u, Vector g) {
		Equation eq = this.getEqnLambda(s_i,a, u, g);
        //Solver solver = new Solver();
        //Vector lmd_solve = solver.solveCGS(eq.A, eq.f);
        
        SolverJBLAS sol = new SolverJBLAS();
		Vector lmd_solve = sol.solveDGESV(eq.A, eq.f);
        return lmd_solve;
    }
	
	protected void setDirichlet(SparseBlockMatrix BM, SparseBlockVector BV,
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
	
	public void imposeDirichletCondition(SparseBlockMatrix BM, SparseBlockVector BV,
			MathFunc diri) {
		ElementList eList = mesh.getElementList();
		//int nNode = mesh.getNodeList().size();
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
						setDirichlet(BM,BV,dof.getGlobalIndex(),diri.apply(v));
						//setDirichlet(BM,BV,nNode+dof.getGlobalIndex(),diri.value(v));
						//setDirichlet(BM,BV,nNode*2+dof.getGlobalIndex(),diri.value(v));
					}
				}
			}
		}
	}
	
	public void imposeDirichletCondition(SparseBlockMatrix BM, SparseBlockVector BV,
			int nDataBlock, MathFunc diri) {
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
							setDirichlet(BM,BV,k*nNode+dof.getGlobalIndex(),diri.apply(v));
							setDirichlet(BM,BV,(nDataBlock+k)*nNode+dof.getGlobalIndex(),diri.apply(v));
						}
						setDirichlet(BM,BV,nDataBlock*2*nNode+dof.getGlobalIndex(),diri.apply(v));
					}
				}
			}
		}
	}
	
	public void imposeDirichletCondition(SparseBlockMatrix BM, SparseBlockVector BV,
			int nDataBlock, MathFunc[] u, MathFunc[] g, MathFunc[] lambda) {
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
					//NodeType.Robin！！！
					if(n.getNodeType() == NodeType.Robin) {
						Variable v = Variable.createFrom(u[0], n, n.globalIndex);
						//循环每次测量
						for(int k=0;k<nDataBlock;k++) {
							setDirichlet(BM,BV,k*nNode+dof.getGlobalIndex(),
									(g[k].apply(v)-u[k].apply(v))/beta
									//0.0
									);
//							setDirichlet(BM,BV,(nDataBlock+k)*nNode+dof.getGlobalIndex(),
//									-lambda[k].value(v)
//									//0.0
//									);
						}
						//setDirichlet(BM,BV,nDataBlock*2*nNode+dof.getGlobalIndex(),0.0);
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
//	public SparseBlockMatrix getSearchDirectionMatrix(Vector ak, Vector uk) {
//		SparseBlockMatrix BM = new SparseBlockMatrix(3,3);
//		Matrix M = this.getM();
//		Matrix A = this.getA(ak,true);
//		Matrix AT = this.getAT(ak);
//		Matrix C = this.getC(ak,uk);
//		Matrix CT = this.getCT(ak,uk);
//		Matrix R = this.getBR();
//		
//		Matrix BM13 = new SparseMatrix(M.getRowDim(),R.getColDim());
//		Matrix BM22 = new SparseMatrix(A.getRowDim(),A.getColDim());
//		Matrix BM31 = new SparseMatrix(R.getRowDim(),M.getColDim());
//		
//		BM.setBlock(1, 1, M);
//		BM.setBlock(1, 2, AT);
//		BM.setBlock(1, 3, BM13);
//		
//		BM.setBlock(2, 1, A);
//		BM.setBlock(2, 2, BM22);
//		BM.setBlock(2, 3, C);
//		
//		BM.setBlock(3, 1, BM31);
//		BM.setBlock(3, 2, CT);
//		BM.setBlock(3, 3, R);
//		
//		return BM;
//	}
	
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
	public SparseBlockMatrix getSearchDirectionMatrix(Vector ak, Vector[] uk,
			Vector[] uk_x, Vector[] uk_y,
			List<ParamOfLightSource> paramList) {
		
		int nDataBlock = uk.length;
		SparseBlockMatrix BM = new SparseBlockMatrix(nDataBlock*2+1,nDataBlock*2+1);
		
		SparseMatrix R = this.getBR();
		BM.setBlock(nDataBlock*2+1, nDataBlock*2+1, R.setName("R"));
		
		for(int i=1;i<=nDataBlock;i++) {
			SparseMatrix M = this.getM();
			SparseMatrix A = this.getA(ak,true);
			SparseMatrix AT = this.getAT(ak);
			//SparseMatrix C = this.getCT(ak,uk[i-1],uk_x[i-1],uk_y[i-1]).trans();
			SparseMatrix C = this.getC(ak,uk[i-1],uk_x[i-1],uk_y[i-1]);
			SparseMatrix CT = this.getCT(ak,uk[i-1],uk_x[i-1],uk_y[i-1]);
			for(int j=1;j<=nDataBlock;j++) {
				if(i==j) {
					BM.setBlock(i, j, M.setName("M"+i));
					BM.setBlock(i, nDataBlock+j, AT.setName("AT"+i));
					BM.setBlock(nDataBlock+i, j, A.setName("A"+i));
					BM.setBlock(nDataBlock+i, nDataBlock*2+1, C.setName("C"+i));
					BM.setBlock(nDataBlock*2+1, nDataBlock+i, CT.setName("CT"+i));
					
					BM.setBlock(nDataBlock+i, nDataBlock+i, 
							new SparseMatrixRowMajor(A.getRowDim(),AT.getColDim()));
					BM.setBlock(i, nDataBlock*2+1, 
							new SparseMatrixRowMajor(M.getRowDim(),R.getColDim()));
					BM.setBlock(nDataBlock*2+1, i, 
							new SparseMatrixRowMajor(R.getRowDim(),M.getColDim()));
					
				} else {
					
					SparseMatrix M0 = new SparseMatrixRowMajor(M.getRowDim(),M.getColDim());
					SparseMatrix AT0 = new SparseMatrixRowMajor(AT.getRowDim(),AT.getColDim());
					SparseMatrix A0 = new SparseMatrixRowMajor(A.getRowDim(),A.getColDim());
					SparseMatrix ATA0 = new SparseMatrixRowMajor(A.getRowDim(),AT.getColDim());
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
	public SparseBlockMatrix changeBlockColumn(SparseBlockMatrix BM, int col1, int col2) {
		int blockRow = BM.getRowBlockDim();
		int blockCol = BM.getColBlockDim();
		
		SparseBlockMatrix newBM = new SparseBlockMatrix(blockRow, blockCol);
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
	
	public SparseBlockMatrix changeBlockColumn(SparseBlockMatrix BM, int nDataBlock) {
		int blockRow = BM.getRowBlockDim();
		int blockCol = BM.getColBlockDim();
		
		SparseBlockMatrix newBM = new SparseBlockMatrix(blockRow, blockCol);
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

	public void testRefine2() {
		ElementList eList = mesh.getElementList();
		ElementList eToRefine = new ElementList();

		for(int i=1;i<=eList.size();i++) {
			Element e = eList.at(i);
			NodeList nodes = e.nodes;
			for(int j=1;j<=nodes.size();j++) {
				Node node = nodes.at(j);
				if(node.coord(1)>2.4 && node.coord(1)<3.8 &&
						node.coord(2)>1.3) {
					eToRefine.add(e);
					break;
				}
			}
		}
		
		this.refineMesh(eToRefine);
	}

	
	/**
	 * 加密网格，会改变mesh和meshBig
	 * @param ak
	 * @param factor
	 */
	public void refineAllMesh(Vector ak, double factor) {
		ElementList eToRefine = Tools.computeRefineElement(mesh, ak, factor);
		writeElementIndex("eToRefineIndex.txt",eToRefine);
		this.refineMesh(eToRefine);
	}
	
	public void refineAllMeshMax(Vector indicator, double factor) {
		ElementList eToRefine = Tools.computeRefineElementMax(mesh, indicator, factor);
		writeElementIndex("eToRefineIndex.txt",eToRefine);
		this.refineMesh(eToRefine);
	}
	
	public void writeElementIndex(String fileName, ElementList eList) {
	    FileOutputStream out = null;
	    PrintWriter br = null;
	    try {
			File file = new File("./"+this.getOutputFolder()+"/"+fileName);
			out = new FileOutputStream(file);
			OutputStreamWriter writer = new OutputStreamWriter(out, "UTF-8");
			br = new PrintWriter(writer);
			for(int i=1;i<=eList.size();i++) {
				br.println(eList.at(i).globalIndex);
			}
			if(br != null)
				br.close();
			if(out != null)
				out.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	public ObjIndex readElementIndex(String fileName) {
		FileInputStream in;
		ObjIndex rlt = new ObjIndex();
		try {
			File file = new File("./"+this.getOutputFolder()+"/"+fileName);
			in = new FileInputStream(file);

			InputStreamReader reader = new InputStreamReader(in,"UTF-8");
			BufferedReader br = new BufferedReader(reader);
	
			String str = null;

			while((str = br.readLine()) != null){
				rlt.add(Integer.parseInt(str.trim()));
			}
			br.close();
			in.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
		return rlt;
	}	
	
	
	/**
	 * Get measurement data of light source s_i
	 * 
	 * @param s_i
	 * @return
	 */
	public Vector getMeasurementData(int s_i) {
		return null;
	}
	
	/**
	 * 先标记网格边界，再置零内部结点值，该方法不会影响原有边界标记
	 * @param mesh
	 * @param v
	 * @return
	 */
	public Vector clearInnerValues(Mesh mesh, Vector v) {
		HashMap<NodeType, MathFunc> mapNTF = new HashMap<NodeType, MathFunc>();
		mapNTF.put(NodeType.Robin, null);
		Map<NodeType, MathFunc> mapNTFOld = mesh.getMarkBorderMap(); 
		
		mesh.clearBorderNodeMark();
		mesh.markBorderNode(mapNTF);
		NodeList nodes = mesh.getNodeList();
		for(int j=1;j<=nodes.size();j++) {
			//System.out.println(nodes.at(j)+"   "+nodes.at(j).getNodeType());
			//Node node = nodes.at(j);
			//if(node.coord(1))
			if(nodes.at(j).getNodeType()==NodeType.Inner)
				v.set(j,0.0);
			
			
		}
		mesh.clearBorderNodeMark();
		mesh.markBorderNode(mapNTFOld);
		
		return v;
	}
	
	public Vector getAFromGCM(Mesh oldMesh, Mesh newMesh) {
		if(newMesh == null) {
			if(this.useVectorMu_a) {//2011/10/18
				return Tools.interplateFrom(aMesh, oldMesh, 
						new Vector2MathFunc(aGuess,aMesh,"x","y").
						setDefaultFunction(FC.c(this.aBackground)));
			}
			return Tools.function2vector(oldMesh, modelGuess.getMu_a());
		} else if(oldMesh != null && newMesh != null){
			//当网格加密后，aGlob采用插值计算出来，而不采用modelGuess.mu_a，
			//否则会导致重构结果在加密单元的中间结点产生小突起（整个看起来有很多小突起）
			if(this.useVectorMu_a) {//2011/10/18
				return Tools.interplateFrom(aMesh, newMesh, 
						new Vector2MathFunc(aGuess,aMesh,"x","y").
						setDefaultFunction(FC.c(this.aBackground)));
			}
			if(aGlob ==  null)
				aGlob = Tools.function2vector(oldMesh, modelGuess.getMu_a());
			return Tools.interplateFrom(oldMesh, newMesh, 
					new Vector2MathFunc(aGlob,oldMesh,"x","y").setDefaultFunction(FC.c(this.aBackground)));
		}
		return null;
	}

	public List<ParamOfLightSource> generateSimulateData() {
		List<ParamOfLightSource> paramList = new ArrayList<ParamOfLightSource>();
		for(int i=0; i<LSx.length; i++) {
			reinitModelLight(i);
			ParamOfLightSource para = new ParamOfLightSource();
			para.s_i = i;
			para.g = solveRealU(i);
			if(debug) {
				meshBig.writeNodesInfo(String.format("%s/meshBig%02d.dat",this.getOutputFolder(),iterNum));
				MathFunc gx = new DuDx(mesh,new Vector2MathFunc(para.g,mesh,"x","y"),"x");
				MathFunc gy = new DuDx(mesh,new Vector2MathFunc(para.g,mesh,"x","y"),"y");
				Tools.plotFunction(mesh, this.getOutputFolder(), 
						String.format("testM%02d_g_grad.dat",i), gx, gy);
			}
			if(!bTestWholdDomain) {
				clearInnerValues(mesh, para.g);
			}
			plotVector(mesh, para.g, String.format("M%02d_g.dat",i));
			paramList.add(para);
		}
		
		//以边界上的测量结果作为边界条件，再以aGlob为系数，求解状态方程，作为“测量”解，在整个区域上计算u-z,
		//而不是只在边界上计算u-g (g=z|_\Gamma)
		if(bTestBoundaryAsWholdDomain) {
			for(int i=0; i<LSx.length; i++) {
				reinitModelLight(i);
				ParamOfLightSource para = paramList.get(i);
				//重新更新para.g
				Vector tmp = para.g;
				para.g = solveStateEquation(aGlob, para.g, null, null);
				plotVector(mesh, para.g, String.format("M%02d_gApproximate.dat",i));
				plotVector(mesh, FMath.axpy(-1.0, tmp, para.g), String.format("M%02d_g_gApproximate_diff.dat",i));
			}
		}
		return paramList;
	}
	
	public List<ParamOfLightSource> getMeasurementData() {
		return null;
	}
	
	
	public Vector gaussNewtonIterateMulti(
			//Input
			List<ParamOfLightSource> paramList, Vector a0,
			//Output
			Vector refineIndicator
			) {
		
		int nDataBlock = paramList.size();
		
		//*************************Initial Values********************
		//u0: 求解状态方程得到
		Vector[] u0 = new SparseVector[nDataBlock];
		Vector[] u0_x = new SparseVector[nDataBlock];
		Vector[] u0_y = new SparseVector[nDataBlock];
		//lambda0: 求解伴随方程得到
		Vector[] lambda0 = new SparseVector[nDataBlock];
		
		for(int i=0;i<nDataBlock;i++) {
			this.reinitModelLight(i);
			ParamOfLightSource param =  paramList.get(i);
			
			//在大区域上求解状态方程，然后计算导数，最后截取到小区域上
			//用于Neumann边界条件情形
			if(this.useVectorMu_a) {//2011/10/18
				Vector aGuessBig = Tools.interplateFrom(aMesh, meshBig, 
						new Vector2MathFunc(aGuess,aMesh,"x","y").
						setDefaultFunction(FC.c(this.aBackground)));
				modelGuess.setMu_a(new Vector2MathFunc(aGuessBig));
			}
			//使用modelGuess计算是不是有问题？应该用mu_a=a0？
			Vector uBig = this.modelGuess.solveNeumann(meshBig);
			plotVector(meshBig, uBig, String.format("M%02d_u0Big.dat",i));
			Vector uBig_x = Tools.computeDerivativeFast(meshBig, uBig, "x");
			Vector uBig_y = Tools.computeDerivativeFast(meshBig, uBig, "y");
			u0_x[i] = Tools.extractData(meshBig, mesh, uBig_x);
			u0_y[i] = Tools.extractData(meshBig, mesh, uBig_y);
			
			//u0初值的边界条件赋予测量数据
			if(bTestWholeDomainDirichletBoundary)
				u0[i] = this.solveStateEquation(a0,param.g, null, null);
			else
				u0[i] = this.solveStateEquation(a0,null, u0_x[i], u0_y[i]);
			plotVector(mesh, u0[i], String.format("M%02d_u0.dat",i));
			plotVector(mesh, FMath.axpy(-1.0, param.g, u0[i]), String.format("M%02d_u_g.dat",i));
			
			
			lambda0[i] = this.solveAdjointEquation(param.s_i, a0, u0[i], param.g);
			mesh.writeNodesInfo(String.format("%s/meshLambda%02d.dat",this.getOutputFolder(),iterNum));

			plotVector(mesh, lambda0[i], String.format("M%02d_lambda0.dat",i));
		}
		//************************************************************
		//**************************TEST***************************
		Matrix stiff = this.testGetCT(a0, u0[0], u0_x[0], u0_y[0]);
		Equation load  = this.testGetCTLoad(u0,u0_x,u0_y,lambda0,a0,FC.C0);
        SolverJBLAS sol2 = new SolverJBLAS();
		Vector testLambda = sol2.solveDGESV(stiff, load.f);
		plotVector(mesh, testLambda, String.format("testLambda.dat"));
		//************************************************************
		
		//迭代求解
		int maxIter = this.maxIterNumPerRefinement[this.refineNum];
		SparseBlockVector f = new SparseBlockVector(nDataBlock*2+1);
		for(int iter=0;iter<maxIter;iter++) {
			
			this.iterNum = iter;
			br.println(">>>>iter="+iter);
			
			//步长矩阵
			SparseBlockMatrix BM = this.getSearchDirectionMatrix(a0, u0, u0_x, u0_y, paramList);
			
			//步长右端向量
			for(int i=1;i<=nDataBlock;i++) {
				this.reinitModelLight(i-1);
				ParamOfLightSource param =  paramList.get(i-1);
				
				//状态方程余量
				SparseVector resLlmd = null;
				if(bTestWholeDomainDirichletBoundary)
					resLlmd = this.getResLlmd(a0, u0[i-1],param.g, null, null).scale(-1.0);
				else
					resLlmd = this.getResLlmd(a0, u0[i-1],null, u0_x[i-1], u0_y[i-1]).scale(-1.0);
		        plotVector(mesh,resLlmd,
		        		String.format("M%02d_Res_Llambda%02d.dat",i-1,this.iterNum));
				f.setBlock(nDataBlock+i, resLlmd);			
				
				//伴随方程余量
				SparseVector resLu = this.getResLu(param.s_i, 
					a0, u0[i-1], param.g,lambda0[i-1]).scale(-1.0);//？1.0
		        plotVector(mesh,resLu,
		        		String.format("M%02d_Res_Lv%02d.dat",i-1,this.iterNum));

				f.setBlock(i, resLu);
			}
			//正则化参数方程余量
			f.setBlock(nDataBlock*2+1, this.getResLq( 
					u0, u0_x, u0_y, lambda0, a0, aGlob,
					FC.c(this.aBackground) //a(x)的边界条件: back ground of a(x)
					).scale(-1.0));

/*			
			//设置边界条件并求解
			//需要交换矩阵BM的第1, 2列，然后设置边界条件，这样使得矩阵A, A'可逆
			//BlockMatrix newBM = this.changeBlockColumn(BM, nDataBlock);
	        //if(this.bTestWholeDomainDirichletBoundary)
			//	this.imposeDirichletCondition(newBM, f, nDataBlock, FC.c0);
			SchurComplementLagrangianSolver solver = 
				new SchurComplementLagrangianSolver(BM, f,mesh);
			BlockVector x = solver.solveMulti();
			//BlockVector x = solver.solve();//Test for one light source
*/
			//this.imposeDirichletCondition(BM, f, FC.c0);
	        if(this.bTestWholeDomainDirichletBoundary)
	        	this.imposeDirichletCondition(BM, f, nDataBlock,FC.C0);
//	        else {
//	        	Function [] fu0 = new Function[nDataBlock];
//	        	Function [] fg = new Function[nDataBlock];
//	        	Function [] flambda0 = new Function[nDataBlock];
//				for(int i=0;i<nDataBlock;i++) {
//					ParamOfLightSource param =  paramList.get(i);
//					fu0[i] = new Vector2Function(u0[i]);
//					flambda0[i] = new Vector2Function(lambda0[i]);
//					fg[i] = new Vector2Function(param.g);
//				}
//	        	this.imposeDirichletCondition(BM, f, nDataBlock,fu0, fg, flambda0);
//	        }
	        
	        //BlockMatrix newBM = this.changeBlockColumn(BM, nDataBlock);
	        System.out.println(String.format("*************Matrix BM(%d by %d)************",BM.getRowDim(),BM.getColDim()));
			SolverJBLAS sol = new SolverJBLAS();
			SparseBlockVector x = sol.solveDGESV(BM, f);

			for(int i=1;i<=nDataBlock;i++) {
				plotVector(mesh, x.getBlock(i), String.format("M%02d_delta_v%02d.dat",i-1,iter));
				plotVector(mesh, x.getBlock(nDataBlock+i), String.format("M%02d_delta_lambda%02d.dat",i-1,iter));
			}
			Vector delta_a = x.getBlock(nDataBlock*2+1);
			plotVector(mesh, delta_a, String.format("delta_a%02d.dat",iter));
			
			//cutX cutY
			//截取一定范围的结果，其他地方都置零（只截取delta_a还不够，还需要重新计算v和lambda!!!)
			NodeList nodes = this.mesh.getNodeList();
			for(int j=1;j<=nodes.size();j++) {
				Node node = nodes.at(j);
				if(node.coord(1)<this.cutX[0] || node.coord(1)>this.cutX[1])
					delta_a.set(j, 0.0);
				else if(node.coord(2)<this.cutY[0] || node.coord(2)>this.cutY[1])
					delta_a.set(j, 0.0);
			}
			plotVector(mesh, delta_a, String.format("delta_cut_a%02d.dat",iter));
			
			//Smooth
			for(int nsm=0;nsm<this.smoothNum;nsm++)
				delta_a = Utils.gaussSmooth(mesh, delta_a, 2, 0.5);
			plotVector(mesh, delta_a, String.format("delta_cut_smooth_a%02d.dat",iter));
			
			//Cut value of delta_a
			double max = FMath.max(delta_a);
			double min = FMath.min(delta_a);
			for(int j=1;j<=nodes.size();j++) {
				double val = delta_a.get(j);
				if(val>=0 && val<this.cutThreshold*max)
					delta_a.set(j, 0.0);
				else if(val<=0 && val>this.cutThreshold*min)
					delta_a.set(j, 0.0);
			}
			plotVector(mesh, delta_a, String.format("delta_cut_smooth_cut_a%02d.dat",iter));

//得到delta_a后，显示求解delta_v
//			for(int i=0;i<nDataBlock;i++) {
//				//!!!
//				this.reinitModelLight(i);
//				
//				ParamOfLightSource param =  paramList.get(i);
//				//????????????????????delta_a
//				Vector dvs = solveDeltaU(a0, delta_a, u0[i]);
//				plotVector(mesh, dvs, String.format("M%02d_rlt_delta_v_solve%02d.dat",i,iterNum));
//			}
			
			//-----------------------步长递减法搜索步长---------------------------
			double stepBase = 0.0;
			double stepDelta = 1.0; //默认初始步长
			
			//如果指定初始步长，则按照指定的开始，否则默认从1.0开始搜索
			if(this.initStepLengthPerRefinement != null &&
					this.initStepLengthPerRefinement.length > this.refineNum)
				stepDelta = this.initStepLengthPerRefinement[this.refineNum];
			
			//计算当前light intensity误差范数：uk-z (如果只是边界： uk|_\Gamma-g)
			double[] norm2_v_g = new double[nDataBlock];
			for(int i=0;i<nDataBlock;i++) {
				ParamOfLightSource param =  paramList.get(i);
				Vector u0Tmp = u0[i];
				if(!bTestWholdDomain) {
					u0Tmp = u0[i].copy();
					clearInnerValues(mesh, u0Tmp);
				}
				norm2_v_g[i] = FMath.axpy(-1.0, param.g, u0Tmp).norm2();
			}
			
			//ak+da后，计算light intensity误差，判断是否减小，以下过程寻找“近似最小”：
			//如果一直减少就往下寻找，直到变大，取上一个步长
			//double lastFindStepLen = -1.0;
			//double lastNormSearch = 0.0;
			for(int searchNum=0;searchNum<maxSearchNum;searchNum++) {
				Vector a0Search = a0.copy();
				a0Search.add(stepBase+stepDelta, delta_a);

				boolean bFind = true;
				br.println(String.format("***Search %02d, stepLength=%f",searchNum,stepBase+stepDelta));
				System.out.println(String.format("Search %02d, stepLength=%f",searchNum,stepBase+stepDelta));
				double[] norm2_v_gSearch = new double[nDataBlock];
				for(int i=0;i<nDataBlock;i++) {
					this.reinitModelLight(i);
					ParamOfLightSource param =  paramList.get(i);
					Vector vsol = null;
					if(bTestWholeDomainDirichletBoundary)
						vsol = this.solveStateEquation(a0Search,param.g,null,null);
					else
						vsol = this.solveStateEquation(a0Search,null,u0_x[i],u0_y[i]);
					plotVector(mesh, vsol, String.format("M%02d_rlt_search_v_solve%02d.dat",i,iterNum));

					Vector vsolTmp = vsol;
					if(!bTestWholdDomain) {
						vsolTmp = vsol.copy();
						clearInnerValues(mesh, vsolTmp);
					}
					norm2_v_gSearch[i] = FMath.axpy(-1.0, param.g, vsolTmp).norm2();
					br.println("\tLS"+i+" NormSearch="+norm2_v_gSearch[i]+"  NormLast="+norm2_v_g[i]);
					System.out.println("LS"+i+" NormSearch="+norm2_v_gSearch[i]+"  NormLast="+norm2_v_g[i]);
					
					
					Vector a1 = FMath.axpy(stepDelta, delta_a, a0);
					double normInf_da = delta_a.normInf();
					double normInf_a0 = a0.normInf();
					double normInf_a1 = a1.normInf();
					if( norm2_v_gSearch[i] > this.maxTargetNormIncFactor*norm2_v_g[i] || 
						stepDelta*normInf_da > normInf_a0 ||
						normInf_a1 > this.maxInfNormIncFactor*normInf_a0) {
						
						if(norm2_v_gSearch[i] > this.maxTargetNormIncFactor*norm2_v_g[i]) 
							br.println(String.format("\tnorm2_v_gSearch[i] > maxTargetNormIncFactor*norm2_v_g[i] (%.5f>%.5f)",norm2_v_gSearch[i],this.maxTargetNormIncFactor*norm2_v_g[i]));
						if(stepDelta*normInf_da>normInf_a0) 
							br.println(String.format("\tstepDelta*normInf_da > normInf_a0 (%.5f>%.5f)",stepDelta*normInf_da,normInf_a0));
						if(normInf_a1 > this.maxInfNormIncFactor*normInf_a0) 
							br.println(String.format("\tnormInf_a1 > this.maxInfNormIncFactor*normInf_a0 (%.5f>%.5f)",normInf_a1,this.maxInfNormIncFactor*normInf_a0));

						br.println("\tgo to next search...\n");
						bFind = false;
						break;
					}
				}

				if(bFind){ // && lastFindStepLen > 0.0 && FMath.max(v_gNormSearch) > lastNormSearch) {
					br.println("found!");
					break;
				}

				//lastNormSearch = FMath.max(v_gNormSearch);
				//lastFindStepLen = stepDelta;
				stepDelta = stepReduceFactor*stepDelta;
				
			}
			stepBase = stepDelta;
			br.println(String.format("stepLength=%f",stepBase));
			br.flush();
			System.out.println(String.format("Iter=%02d==========================>stepLength=%f",iter,stepBase));
			
//--------------按比例搜索步长---------------------------			
//			for(int i=1;i<=nDataBlock;i++) {
//				ParamOfLightSource param =  paramList.get(i-1);
//				Vector uk = u0[i-1];
//				if(!bTestWholdDomain) {
//					uk = u0[i-1].copy();
//					clearInnerValues(mesh, uk);
//				}
//				Vector v_z = FMath.axpy(-1.0, param.g, uk);
//				errorNorm = v_z.norm2();
//				br.format("M%02d Error norm2(u-g)=%f \n", i, errorNorm);
//				Vector delta_v = x.getBlock(i);
//				if(!bTestWholdDomain) {
//					delta_v = x.getBlock(i).copy();
//					clearInnerValues(mesh,delta_v);
//				}
//				//stepLength*||delta_v|| < ||v_z||
//				double tmp = errorNorm/delta_v.norm2();
//				if(stepLength > tmp)
//					stepLength = tmp;
//			}
//			if(stepLength * delta_a.normInf() > a0.normInf()*0.2) {
//				stepLength = (0.5/Math.pow(2, iter))*a0.normInf()/delta_a.normInf();
//				//stepLength = (0.1/(iter+1))*a0.normInf()/delta_a.normInf();
//				System.out.println("Iter"+iter+"++++++++++++++++++++++++++++++++=>stepLength="+stepLength);
//				br.println("stepLength_m2="+stepLength);
//			} else {
//				System.out.println("Iter"+iter+"=================================>stepLength="+stepLength);
//				br.println("stepLength_m1="+stepLength);
//			}
//			br.println();
//			if(stepLength > lastStepLength) {
//				stepLength = 0.1*lastStepLength;
//			}
//			lastStepLength = stepLength;
//			errorNorm = 2.0;
			
			
			//这一步更新u0,lambda0实际没有使用，后面是根据更新的a0计算新的u0,lambda0
			for(int i=1;i<=nDataBlock;i++) {
				u0[i-1].add(stepDelta*beta, x.getBlock(i));
				//注意：-stepLength应该是正的
				lambda0[i-1].add(stepDelta*beta, x.getBlock(nDataBlock+i));
			}
			a0.add(stepDelta, delta_a);
			
			for(int i=1;i<=nDataBlock;i++) {
				plotVector(mesh, u0[i-1], String.format("M%02d_rlt_v%02d.dat",i-1,iter));
				plotVector(mesh, lambda0[i-1], String.format("M%02d_rlt_lambda%02d.dat",i-1,iter));
			}
			plotVector(mesh, a0, String.format("rlt_a%02d_refine%02d.dat",iter,refineNum));
			
			//计算新a0对应的v，再计算lmd，用来验证
			for(int i=0;i<nDataBlock;i++) {
				this.reinitModelLight(i);
				ParamOfLightSource param =  paramList.get(i);
				
				Vector vsol = null;
				if(bTestWholeDomainDirichletBoundary)
					vsol = this.solveStateEquation(a0,param.g,null,null);
				else
					vsol = this.solveStateEquation(a0,null,u0_x[i],u0_y[i]);
				plotVector(mesh, vsol, String.format("M%02d_rlt_v_solve%02d.dat",i,iterNum));
				plotVector(mesh, FMath.axpy(-1.0, vsol, u0[i]), String.format("M%02d_rlt_v_diff%02d.dat",i,iterNum));
				//直接计算u0+du
				u0[i] = vsol;
				
				Vector lamsol = this.solveAdjointEquation(param.s_i, a0, vsol, param.g);
				plotVector(mesh, lamsol, String.format("M%02d_rlt_lambda_solve%02d.dat",i,iterNum));
				plotVector(mesh, FMath.axpy(-1.0, lamsol, lambda0[i]), String.format("M%02d_rlt_lambda_diff%02d.dat",i,iterNum));
				//直接计算lambda0+dl
				lambda0[i] = lamsol;
				
			}
			
			//DWR refine indicator
			Vector[] laplaceLambda0 = new SparseVector[nDataBlock];
			Vector[] laplaceU0 = new SparseVector[nDataBlock];
			Vector[] multLU = new SparseVector[nDataBlock];
			for(int i=0;i<nDataBlock;i++) {
				laplaceLambda0[i] = Tools.computeLaplace2D(mesh, lambda0[i]);
				laplaceU0[i] = Tools.computeLaplace2D(mesh, u0[i]);
				multLU[i] = FMath.axMuly(1.0, laplaceLambda0[i], laplaceU0[i]);
				plotVector(mesh,laplaceLambda0[i],String.format("M%02dlaplaceLambda_%02d.dat",i,iterNum));
				plotVector(mesh,laplaceU0[i],String.format("M%02dlaplaceU0_%02d.dat",i,iterNum));
			}
			Vector ind = FMath.sum(multLU);
			plotVector(mesh,ind,String.format("DWR_indicator_%02d.dat",iterNum));
			refineIndicator.set(FMath.abs(ind));
			
			//if(errorNorm < 1.0) break;
			
		}
		return a0;

		
	}	
	
	protected void work(Vector a0, int beginRefineNum, boolean simulate) {
		beginLog();
		//Begin refinement iteration
		Vector ak = a0;
		if(simulate) {
			for(int i=beginRefineNum; i<totalRefineNum; i++) {
				refineNum = i;
				
				//Write mesh file in local folder for later use in 'reStart()' function if necessary
				MeshWriter.write2DMesh(mesh, String.format("./%s/%s",getOutputFolder(),gridFileSmall));
				MeshWriter.write2DMesh(meshBig, String.format("./%s/%s",getOutputFolder(),gridFileBig));
				
				List<ParamOfLightSource> paramList = generateSimulateData();
				Vector refineIndicator = new SparseVectorHashMap();
				
				//Gauss-Newtom Iteration
				Vector aNew = gaussNewtonIterateMulti(paramList,ak,refineIndicator);
				
				//Refine all meshes (mesh and meshBig)
				Mesh oldMesh = mesh.copy();
				//refineAllMesh(aNew,refineFactors[i]);
				refineAllMeshMax(refineIndicator,refineFactorsForMesh[i]);
				//mesh.printMeshInfo();
				mesh.writeNodesInfo(getOutputFolder()+"/meshInfo.dat");
				
				//Interplate ak from old mesh to new refined mesh
				ak = Tools.interplateFrom(oldMesh,mesh,
						new Vector2MathFunc(aNew,oldMesh,"x","y").setDefaultFunction(FC.c(this.aBackground)));
				//Interplate aGlob from old mesh to new refined mesh
				//or
				//Read a(x) from GCM (Global Convergence Method) based on new refined mesh
				
				//2011-10-24 设置aGlob为ak
				this.aGlob = getAFromGCM(oldMesh, this.mesh);
				//this.aGlob = ak.copy();
				
				
				//Set new output folder index
				setOutputFolderIndex(i+1);
				
				String oFolder = getOutputFolder()+"/Input";
			    if(!oFolder.isEmpty()) {
				    File file = new File(oFolder);
					if(!file.exists()) {
						file.mkdirs();
					}
			    }
			    
				//Plot parameters: a0, aReal, aGlob (After refinement)
				plotVector(mesh, ak, "a0.dat");
				if(this.useVectorMu_a) {//2011/10/18
					//mesh already been refined now!
					Vector aRealRefine = Tools.interplateFrom(aMesh, mesh, 
							new Vector2MathFunc(aReal,aMesh,"x","y").
							setDefaultFunction(FC.c(this.aBackground)));
					modelReal.setMu_a(new Vector2MathFunc(aRealRefine));
				}
				plotFunction(mesh, modelReal.getMu_a(), String.format("aReal_refine%02d.dat",i));
				if(this.useVectorMu_a) {//2011/10/18
					plotFunction(mesh, modelReal.getMu_a(), String.format("Input/input_real_mu_a.dat",i));
					//mesh already been refined now!
					//重新运行已经不需要aGuess
					Vector aRealRefine = Tools.interplateFrom(aMesh, mesh, 
							new Vector2MathFunc(aGuess,aMesh,"x","y").
							setDefaultFunction(FC.c(this.aBackground)));
					modelGuess.setMu_a(new Vector2MathFunc(aRealRefine));
					plotFunction(mesh, modelGuess.getMu_a(), String.format("Input/input_guess_mu_a.dat",i));
				}

				plotVector(mesh, aGlob, "aGlob.dat");
				plotVector(mesh, FMath.axpy(-1.0, ak, 
						Tools.function2vector(mesh, modelReal.getMu_a())),"aReal_diff.dat");			
				
				if(this.useVectorMu_a) {//2011/10/18
					//mesh already be refined now!
					Vector aRealBigRefine = Tools.interplateFrom(aMesh, meshBig, 
							new Vector2MathFunc(aReal,aMesh,"x","y").
							setDefaultFunction(FC.c(this.aBackground)));
					modelReal.setMu_a(new Vector2MathFunc(aRealBigRefine));
				}
				plotFunction(meshBig, modelReal.getMu_a(), String.format("aRealBig_refine%02d.dat",i));
			}
		} else {
			//TODO
		}
		endLog();
		
	}
	
	/**
	 * 从中间结果回复计算，以上次加密编号和上次计算迭代次数的重构结果开始新的一轮加密计算
	 * 网格不从根目录读取，而是从产生的子目录读取
	 * 
	 * @param lastRefineNum 上次加密编号
	 * @param lastIterationNum 上次计算迭代次数
	 * @param simulate
	 */
	public void reStart(int lastRefineNum, int lastIterationNum, boolean simulate) {
		//Update grid files to refined files
		setOutputFolderIndex(0);
		readMesh("./"+getOutputFolder());
		for(int i=0;i<lastRefineNum;i++) {
			ObjIndex eleIndex = this.readElementIndex("eToRefineIndex.txt");
			ElementList eToRefine = new ElementList();
			for(int j=0;j<eleIndex.size();j++) {
				eToRefine.add(this.mesh.getElementList().at(eleIndex.get(j)));
			}
			this.refineMesh(eToRefine);
			this.mesh.writeNodesInfo(String.format("_%d.dat", i));
			this.meshBig.writeNodesInfo(String.format("_b%d.dat", i));
			setOutputFolderIndex(i+1);
		}
		
		//读取配置有问题：如果系数是向量，需要读入一个网格，这个网格何时读取是问题
		init();

		setOutputFolderIndex(lastRefineNum);
		//Read initial guess of a0 （上次的计算结果）
		Vector a0 = DataReader.readVector(String.format("./%s/rlt_a%02d_refine%02d.dat",
				getOutputFolder(),lastIterationNum,lastRefineNum));
		
		//Refine all meshes (mesh and meshBig)
		Mesh oldMesh = mesh.copy();
		refineAllMesh(a0,refineFactorsForMesh[lastRefineNum]);
		//mesh.printMeshInfo();
//        HashMap<NodeType, Function> mapNTF =
//            new HashMap<NodeType, Function>();
//	    mapNTF.put(NodeType.Dirichlet, null);
//	    mesh.markBorderNode(mapNTF);
        //mesh.writeNodesInfo(String.format("mesh%02d.dat",iterNum));

		//Interplate ak from old mesh to new refined mesh
		a0 = Tools.interplateFrom(oldMesh,mesh,
				new Vector2MathFunc(a0,oldMesh,"x","y").setDefaultFunction(FC.c(this.aBackground)));
		//Interplate aGlob from old mesh to new refined mesh
		//or
		//Read a(x) from GCM (Global Convergence Method) based on new mesh
		if(simulate) {
			this.aGlob = getAFromGCM(oldMesh,this.mesh);
		} else {
			//TODO
		}
		
		//Set new output folder index
		setOutputFolderIndex(lastRefineNum+1);
		
		//Plot parameters: a0, aReal, aGlob (After refinement)
		plotVector(mesh, a0, "a0.dat");
		if(this.useVectorMu_a) {//2011/10/18
			//mesh already be refined now!
			Vector aRealRefine = Tools.interplateFrom(aMesh, mesh, 
					new Vector2MathFunc(aReal,aMesh,"x","y").
					setDefaultFunction(FC.c(this.aBackground)));
			modelReal.setMu_a(new Vector2MathFunc(aRealRefine));
		}
		plotFunction(mesh, modelReal.getMu_a(), String.format("aReal_refine%02d.dat",lastRefineNum+1));
		plotVector(mesh, aGlob, "aGlob.dat");
		plotVector(mesh, FMath.axpy(-1.0, a0, 
				Tools.function2vector(mesh, modelReal.getMu_a())),"aReal_diff.dat");
		
		if(this.useVectorMu_a) {//2011/10/18
			//mesh already be refined now!
			Vector aRealBigRefine = Tools.interplateFrom(aMesh, meshBig, 
					new Vector2MathFunc(aReal,aMesh,"x","y").
					setDefaultFunction(FC.c(this.aBackground)));
			modelReal.setMu_a(new Vector2MathFunc(aRealBigRefine));
		}
		plotFunction(meshBig, modelReal.getMu_a(), String.format("aRealBig_refine%02d.dat",lastRefineNum+1));
		
		work(a0,lastRefineNum+1,simulate);
	}
	
	/**
	 * 开始计算，网格从根目录读取
	 * 
	 * @param simulate
	 */
	public void start(boolean simulate) {
		readMesh("."); //Read mesh in current path
		//testRefine();
		//testRefine2();
		//mesh.printMeshInfo();
		
		init();
		
		//Read a(x) from GCM (Global Convergence Method)
		this.aGlob = this.getAFromGCM(this.mesh, null);
		
		//Choose initial guess of a0
		Vector a0 = null;
		//a0 = aGlob.copy();
		if(this.useVectorMu_a) {//2011/10/18
			a0 = this.aGuess.copy();
			//modelInit.setMu_a(new Vector2Function(aInit));
			//a0 = Tools.function2vector(mesh, modelInit.getMu_a());
		} else {
			a0 = Tools.function2vector(mesh, modelInit.getMu_a());
		}
		//Plot parameters: a0, aReal, aGlob
		plotVector(mesh, a0, "a0.dat");
		if(this.useVectorMu_a)//2011/10/18
			modelReal.setMu_a(new Vector2MathFunc(aReal));
		plotFunction(mesh,   modelReal.getMu_a(),  String.format("aReal.dat"));
		plotVector(mesh, aGlob, "aGlob.dat");
		plotVector(mesh, FMath.axpy(-1.0, a0, 
				Tools.function2vector(mesh, modelReal.getMu_a())),"aReal_diff.dat");
		
		if(this.useVectorMu_a) {//2011/10/18
			Vector aRealBig = Tools.interplateFrom(aMesh, meshBig, 
					new Vector2MathFunc(aReal,aMesh,"x","y").
					setDefaultFunction(FC.c(this.aBackground)));
			modelReal.setMu_a(new Vector2MathFunc(aRealBig));
		}
		plotFunction(meshBig, modelReal.getMu_a(),String.format("aRealBig.dat"));
		
		work(a0,0,simulate);
	}
	
	
	/**
	 * 默认配置文件名称为VGNDOT.conf，可以通过命令行参数指定配置文件名称。
	 * 结算结果输出结果到子文件夹中，文件夹名称在配置文件中设置。
	 * 输出的子文件可以有多个，按照加密层次编号，从初始网格00开始编号。
	 * 
	 * @param args
	 */
	public static void main(String[] args) {
		System.out.println("Java version: "+System.getProperty("java.version"));
		System.out.println("Usage: args[0]=configFileName (default=VGNDOT.conf)");
		
		VariationGaussNewtonDOTGeneral vgn = new VariationGaussNewtonDOTGeneral();
		
		String configFile = "VGNDOT.conf";
		if(args.length == 1)
			configFile = args[0];
		
		//Read parameters from configFile
		System.out.println("Config file:"+configFile);
		vgn.readParameters(configFile);
		
		vgn.printParameters();
		
		//2011-10-24 设置aGlob为ak
		//2011-10-26 取消设置aGlob为ak
		vgn.start(true);
		
		//vgn.reStart(1,0,true);
	}
}
