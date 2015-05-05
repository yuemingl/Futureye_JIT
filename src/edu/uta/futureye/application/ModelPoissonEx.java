package edu.uta.futureye.application;

import java.util.HashMap;

import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.algebra.solver.external.SolverJBLAS;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.NodeType;
import edu.uta.futureye.core.intf.Assembler;
import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.io.MeshReader;
import edu.uta.futureye.lib.assembler.AssemblerScalar;
import edu.uta.futureye.lib.element.FEBilinearRectangle;
import edu.uta.futureye.lib.weakform.WeakFormLaplace2D;
import edu.uta.futureye.util.container.ElementList;

/**
 * Solve a Poisson's equation (k=1,c=0)
 * 
 * Model: -k*\Delta{u} + c*u = f
 * 
 * @author liuyueming
 *
 */
public class ModelPoissonEx {
	
	public Mesh mesh;
	
	//Right hand side
	public MathFunc f = null;
	
	//Coefficient
	public MathFunc k = FC.C1;
	public MathFunc c = FC.C0;
	
	public ModelPoissonEx(Mesh mesh) {
		this.mesh = mesh;
	}

	/**
	 * 求解混合问题，需要提供函数diriBoundaryMark来标记Dirichlet边界类型，
	 * 其余边界为Neumann类型。
	 *   如果是纯Neumann:
	 *     diriBoundaryMark=null
	 *     diri=null
	 *   如果是纯Dirichlet:
	 *     diriBoundaryMark=null
	 *     diri!=null
	 *
	 * 注意：该方法会改变mesh的边界类型，当使用同一个mesh对象求解多个不同边界条件的问题时，
	 *      需要特别设置对应的边界类型
	 * 
	 * @param mesh
	 * @param diriBoundaryMark: the function that marks which segment on the boundary is Dirichlet boundary 
	 * @param diri: the values of Dirichlet condition
	 * @return
	 */
	public Vector solveMixedBorder(MathFunc diriBoundaryMark, MathFunc diri,
			MathFunc robin_g, MathFunc robin_d) {
		//Mark border type
		HashMap<NodeType, MathFunc> mapNTF = new HashMap<NodeType, MathFunc>();
		if(diriBoundaryMark == null && diri == null) {
			mapNTF.put(NodeType.Robin, null);
		} else if(diriBoundaryMark == null && diri != null) {
			mapNTF.put(NodeType.Dirichlet, null);
		} else {
			mapNTF.put(NodeType.Dirichlet, diriBoundaryMark);
			mapNTF.put(NodeType.Robin, null);
		}
		mesh.clearBorderNodeMark();
		mesh.markBorderNode(mapNTF);
		
		WeakFormLaplace2D weakForm = new WeakFormLaplace2D();
		
		//Right hand side
		weakForm.setF(this.f);


		//Model: \Delta{u} = f
		weakForm.setParam(
				k, c, robin_g, robin_d 
			);
		
		//bugfix 2011-5-7两种方式结果不一样？
		//Assembler assembler = new AssemblerScalarFast(mesh, weakForm);
		Assembler assembler = new AssemblerScalar(mesh, weakForm);
		System.out.println("Begin Assemble...solveMixedBorder");
		assembler.assemble();
		Matrix stiff = assembler.getStiffnessMatrix();
		Vector load = assembler.getLoadVector();
		//Dirichlet condition
		if(diri != null)
			assembler.imposeDirichletCondition(diri);
		System.out.println("Assemble done!");

		//Solver solver = new Solver();
		//Vector u = solver.solveCGS(stiff, load);
        SolverJBLAS sol = new SolverJBLAS();
		Vector u = sol.solveDGESV(stiff, load);
		
		return u;
	}
	
	public Vector solveNeumann() {
		return solveMixedBorder(null,null,null,FC.C1);
	}

	public Vector solveDirichlet(MathFunc diri) {
		return solveMixedBorder(null,diri,null,null);
	}	
	
	public static void gen1(int type) {
		String outputFolder = "ModelPoissonEx";
		//String gridFileSmall = "prostate_test13.grd";
		String gridFileSmall = "prostate_test14.grd";

		MeshReader readerGCM = new MeshReader(gridFileSmall);
		Mesh meshSmall = readerGCM.read2DMesh();
		
		ModelPoissonEx model = new ModelPoissonEx(meshSmall);
		if(type==1)
			model.f = VariationGaussNewtonDOTGeneral.generateTestRealMu_a(0.4, 0.1).M(15);
		else if(type==2)
			model.f = VariationGaussNewtonDOTGeneral.generateTestGuessMu_a(0.4, 0.1).M(15);
		model.c = FC.c(10.0);
		
		
		//Use element library to assign degree of freedom (DOF) to element
		FEBilinearRectangle fe = new FEBilinearRectangle();
		ElementList eList = meshSmall.getElementList();
		for(int i=1;i<=eList.size();i++)
			fe.assignTo(eList.at(i));
		meshSmall.computeNodeBelongsToElements();
		meshSmall.computeNeighborNodes();		
		
		Vector u = model.solveDirichlet(FC.c(0.1));
		if(type==1)
			Tools.plotVector(meshSmall, outputFolder, "uReal.dat", u);
		else if(type==2)
			Tools.plotVector(meshSmall, outputFolder, "uGuess.dat", u);
	}
	
	public static void gen2(int type) {
		String outputFolder = "ModelPoissonEx";
		//String gridFileSmall = "prostate_test13.grd";
		String gridFileSmall = "prostate_test14.grd";

		MeshReader readerGCM = new MeshReader(gridFileSmall);
		Mesh meshSmall = readerGCM.read2DMesh();
		
		ModelPoissonEx model = new ModelPoissonEx(meshSmall);
		model.k = FC.c(0.07);
		if(type==1)
			model.f = VariationGaussNewtonDOTGeneral.generateTestRealMu_a2(meshSmall, 0.1).M(10);
		else if(type==2)
			model.f = VariationGaussNewtonDOTGeneral.generateTestGuessMu_a2(meshSmall, 0.1).M(10);
		model.c = FC.c(10.0);
		
		//Use element library to assign degree of freedom (DOF) to element
		FEBilinearRectangle fe = new FEBilinearRectangle();
		ElementList eList = meshSmall.getElementList();
		for(int i=1;i<=eList.size();i++)
			fe.assignTo(eList.at(i));
		meshSmall.computeNodeBelongsToElements();
		meshSmall.computeNeighborNodes();		
		
		Vector u = model.solveDirichlet(FC.c(0.1));
		if(type==1)
			Tools.plotVector(meshSmall, outputFolder, "uReal2.dat", u);
		else if(type==2)
			Tools.plotVector(meshSmall, outputFolder, "uGuess2.dat", u);
	}
	
	public static void genTest15(int type) {
		String outputFolder = "ModelPoissonEx";
		//String gridFileSmall = "prostate_test13.grd";
		String gridFileSmall = "prostate_test15.grd";

		MeshReader readerGCM = new MeshReader(gridFileSmall);
		Mesh meshSmall = readerGCM.read2DMesh();
		
		ModelPoissonEx model = new ModelPoissonEx(meshSmall);
		model.k = FC.c(0.2);
		double scaleFactor = 40.0;//20.0;//v1-3=8.5,10
		if(type==1)
			model.f = VariationGaussNewtonDOTGeneral.generateRealMu_aTest15(meshSmall, 1/scaleFactor).M(scaleFactor);
		else if(type==2)
			model.f = VariationGaussNewtonDOTGeneral.generateGuessMu_aTest15(meshSmall, 1/scaleFactor).M(scaleFactor);
		model.c = FC.c(10.0);
		
		//Use element library to assign degree of freedom (DOF) to element
		FEBilinearRectangle fe = new FEBilinearRectangle();
		ElementList eList = meshSmall.getElementList();
		for(int i=1;i<=eList.size();i++)
			fe.assignTo(eList.at(i));
		meshSmall.computeNodeBelongsToElements();
		meshSmall.computeNeighborNodes();		
		
		Vector u = model.solveDirichlet(FC.c(0.1));
//		for(int i=1;i<=u.getDim();i++) {
//			u.set(i,u.get(i)+0.005*Math.random());
//		}
		if(type==1)
			Tools.plotVector(meshSmall, outputFolder, "uReal_Test15.dat", u);
		else if(type==2)
			Tools.plotVector(meshSmall, outputFolder, "uGuess_Test15.dat", u);
	}
	
	public static void genTest15Direct() {
		String outputFolder = "ModelPoissonEx";
		String gridFileSmall = "prostate_test15.grd";

		MeshReader readerGCM = new MeshReader(gridFileSmall);
		Mesh meshSmall = readerGCM.read2DMesh();
		MathFunc real = VariationGaussNewtonDOTGeneral.generateRealMu_aTest15(meshSmall, 0.1);
		MathFunc guess = VariationGaussNewtonDOTGeneral.generateGuessMu_aTest15(meshSmall, 0.1);
		Tools.plotFunction(meshSmall, outputFolder, "uReal_Test15.dat", real);
		Tools.plotFunction(meshSmall, outputFolder, "uGuess_Test15.dat", guess);
		
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		//gen1(1);
		//gen2(2);
		//genTest15(1);
		//genTest15(2);
		genTest15Direct();
	}		
}
