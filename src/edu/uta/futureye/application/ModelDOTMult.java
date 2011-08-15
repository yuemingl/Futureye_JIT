package edu.uta.futureye.application;

import java.util.HashMap;

import edu.uta.futureye.algebra.Solver;
import edu.uta.futureye.algebra.SolverJBLAS;
import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.NodeType;
import edu.uta.futureye.core.intf.Assembler;
import edu.uta.futureye.function.AbstractFunction;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.DuDn;
import edu.uta.futureye.function.basic.DuDx;
import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.basic.FDelta;
import edu.uta.futureye.function.basic.Vector2Function;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.function.operator.FMath;
import edu.uta.futureye.io.MeshReader;
import edu.uta.futureye.lib.assembler.AssemblerScalar;
import edu.uta.futureye.lib.element.FEBilinearRectangle;
import edu.uta.futureye.lib.element.FELinearTriangle;
import edu.uta.futureye.lib.weakform.WeakFormLaplace2D;
import edu.uta.futureye.util.Constant;
import edu.uta.futureye.util.container.ElementList;
import edu.uta.futureye.util.container.NodeList;

/**
 * Solver the following model problem:
 * 
 *   -\nabla{(1/(a*k))*\nabla{u}} + u = \delta/a,  in \Omega_0
 *   
 *   -\nabla{(1/(a*k))*\nabla{u}} + u = 0          in \Omega
 * 
 * where 
 *   k = 3*mu_s'
 *   a = a(x) = mu_a(x)
 *   
 * @author liuyueming
 *
 */
public class ModelDOTMult {
	//Light source
	public Function delta = null;
	public Variable lightPosition = null; //light source position
	public int lightNum = -1;
	
	//Inclusion mu_a
	public Function mu_a = null;
	
	//mu_s'
	public Function mu_s = new FC(50.0/3.0);
	
	/**
	 * type=1: one inclusion
	 * type=2: two inclusion
	 * @param incX
	 * @param incY
	 * @param incR
	 * @param maxMu_a
	 * @param type
	 */
	public void setMu_a(double incX, double incY, double incR, double maxMu_a,
			int type) {
		final double fcx = incX;
		final double fcy = incY;
		final double fcr = incR;
		final double fmu_a = maxMu_a;
		final double distance = 0.8;
		if(type == 1) {
			mu_a = new AbstractFunction("x","y"){
				@Override
				public double value(Variable v) {
					double bk = 0.1;
					double dx = v.get("x")-fcx;
					double dy = v.get("y")-fcy;
					if(Math.sqrt(dx*dx+dy*dy) < fcr) {
						double r = fmu_a*Math.cos((Math.PI/2)*Math.sqrt(dx*dx+dy*dy)/fcr); 
						return r<bk?bk:r;
					}
					else
						return bk;
				}
			};
		} else if(type == 2) {
			mu_a = new AbstractFunction("x","y"){
				@Override
				public double value(Variable v) {
					double bk = 0.1;
					double dx = v.get("x")-fcx;
					double dy = v.get("y")-fcy;
					double dx1 = v.get("x")-(fcx+distance);
					double dy1 = v.get("y")-fcy;
					if(Math.sqrt(dx*dx+dy*dy) < fcr) {
						double r = fmu_a*Math.cos((Math.PI/2)*Math.sqrt(dx*dx+dy*dy)/fcr); 
						return r<bk?bk:r;
					}
					else if(Math.sqrt(dx1*dx1+dy1*dy1) < fcr) {
						double r = fmu_a*Math.cos((Math.PI/2)*Math.sqrt(dx1*dx1+dy1*dy1)/fcr); 
						return r<bk?bk:r;
					}
					else
						return bk;
				}
			};			
		}
	}
	
	public void setDelta(double x,double y) {
		this.lightPosition = new Variable();
		this.lightPosition.set("x", x);
		this.lightPosition.set("y", y);
		delta = new FDelta(this.lightPosition,0.01,2e5);
		delta = delta.D(this.mu_a);
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
	public Vector solveMixedBorder(Mesh mesh, 
			Function diriBoundaryMark, Function diri,
			Function robinQ, Function robinD) {
		//Mark border type
		HashMap<NodeType, Function> mapNTF = new HashMap<NodeType, Function>();
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
		weakForm.setF(this.delta);


		//Model: \nabla{1/(3*mu_s'*mu_a)*\nabla{u}} + u = \delta/mu_a
		weakForm.setParam(
				FC.c1.D(mu_s.M(mu_a).M(3.0)), 
				FC.c1, 
				robinQ, 
				robinD //FC.c1.D(mu_s.M(mu_a).M(3.0) : d==k,q=0 (即：u_n + u =0)
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
		//Tools.plotVector(mesh,"",String.format("x.dat"),x);
		
		return u;
	}
	
	public Vector solveNeumann(Mesh mesh) {
		return solveMixedBorder(mesh,null,null,null,FC.c1.D(mu_s.M(mu_a).M(3.0)));
	}

	public Vector solveDirichlet(Mesh mesh, Function diri) {
		return solveMixedBorder(mesh,null,diri,null,null);
	}	
	
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		String outputFolder = "ModelDOTMult";
//		String gridFileBig = "prostate_test3_ex.grd";
//		String gridFileSmall = "prostate_test3.grd";
		String gridFileBig = "prostate_test7_ex.grd";
		String gridFileSmall = "prostate_test7.grd";
//		String gridFileBig = "prostate_test8_ex.grd";
//		String gridFileSmall = "prostate_test8.grd";


		ModelDOTMult model = new ModelDOTMult();
//		model.setMu_a(2.0, 2.5, 0.5, //(x,y;r)
//				0.4, //maxMu_a
//				1); //type
		model.setMu_a(3.0, 2.30, 0.6,
				0.8, //peak value of mu_a
				1); //Number of inclusions
		model.setDelta(2.5, 3.5);
	
		MeshReader readerForward = new MeshReader(gridFileBig);
		Mesh meshBig = readerForward.read2DMesh();
		MeshReader readerGCM = new MeshReader(gridFileSmall);
		Mesh meshSmall = readerGCM.read2DMesh();
		
		Tools.plotFunction(meshBig, outputFolder, "aReal.dat", model.mu_a);
		//Use element library to assign degree of freedom (DOF) to element
		ElementList eList = meshBig.getElementList();
		//FELinearTriangle fe = new FELinearTriangle();
		FEBilinearRectangle fe = new FEBilinearRectangle();
		for(int i=1;i<=eList.size();i++)
			fe.assignTo(eList.at(i));
		meshBig.computeNodeBelongsToElements();
		meshBig.computeNeighborNodes();
		
		eList = meshSmall.getElementList();
		for(int i=1;i<=eList.size();i++)
			fe.assignTo(eList.at(i));
		meshSmall.computeNodeBelongsToElements();
		meshSmall.computeNeighborNodes();		
		
		
		
		//TEST 1.
		Vector uBig = model.solveNeumann(meshBig);
		Tools.plotVector(meshBig, outputFolder, "u_big.dat", uBig);

		//TEST 2.
		Vector uSmallExtract = Tools.extractData(meshBig, meshSmall, uBig);
		Tools.plotVector(meshSmall, outputFolder, "u_small_extract.dat", uSmallExtract);
		Vector uSmallDiri = model.solveDirichlet(meshSmall, new Vector2Function(uSmallExtract));
		Tools.plotVector(meshSmall, outputFolder, "u_small_diri.dat", uSmallDiri);
		Tools.plotVector(meshSmall, outputFolder, "u_small_extract_diri_diff.dat", 
				FMath.axpy(-1.0, uSmallDiri, uSmallExtract));
		
		
		Vector uSmallDiriBoundary = uSmallDiri.copy();
		NodeList nodes = meshSmall.getNodeList();
		for(int i=1;i<=nodes.size();i++) {
			if(nodes.at(i).isInnerNode())
				uSmallDiriBoundary.set(i,0.0);
		}
		
		//不是很准确
//		Vector uSmall_x = Tools.computeDerivativeFast(meshSmall, uSmallExtract, "x");
//		Vector uSmall_y = Tools.computeDerivativeFast(meshSmall, uSmallExtract, "y");
		//整个区域要准确很多
		Vector uBig_x = Tools.computeDerivativeFast(meshBig, uBig, "x");
		Vector uBig_y = Tools.computeDerivativeFast(meshBig, uBig, "y");
		Vector uSmall_x = Tools.extractData(meshBig, meshSmall, uBig_x);
		Vector uSmall_y = Tools.extractData(meshBig, meshSmall, uBig_y);
		DuDn dudn = new DuDn(
				new Vector2Function(uSmall_x),
				new Vector2Function(uSmall_y),
				null);		
		//单元上分片常数=小区域直接求导数
//		Vector2Function fuSmall = new Vector2Function(uSmallExtract,meshSmall,"x","y");
//		DuDx uSmall_x = new DuDx(meshSmall,fuSmall,"x");
//		DuDx uSmall_y = new DuDx(meshSmall,fuSmall,"y");
//		DuDn dudn = new DuDn(
//				uSmall_x,
//				uSmall_y,
//				null);	
		
		Vector uSmallRobin = model.solveMixedBorder(meshSmall, 
				null, null, 
				//new Vector2Function(uSmallDiriBoundary).M(0.5), null);
				dudn.M(FC.c1.D(model.mu_s.M(model.mu_a).M(3.0))), 
				null);//FC.c1.D(model.mu_s.M(model.mu_a).M(3.0)));
		Tools.plotVector(meshSmall, outputFolder, "u_small_robin.dat", uSmallRobin);
		Tools.plotVector(meshSmall, outputFolder, "u_small_extract_robin_diff.dat", 
				FMath.axpy(-1.0, uSmallRobin, uSmallExtract));
		
		
		//First Guess a(x)
//		model.setMu_a(2.2, 2.5, 0.5, //(x,y;r)
//				0.4, //maxMu_a
//				1); //type
		model.setMu_a(3.2, 2.10, 0.6,
				0.8, //peak value of mu_a
				1); //Number of inclusions
		
		Tools.plotFunction(meshBig, outputFolder, "aRealGuess.dat", model.mu_a);
		Vector uBigGuess = model.solveNeumann(meshBig);
		Tools.plotVector(meshBig, outputFolder, "u_big_guess.dat", uBigGuess);
		
		Vector uSmallGuess = Tools.extractData(meshBig, meshSmall, uBigGuess);
		Tools.plotVector(meshSmall, outputFolder, "u_small_guess.dat", uSmallGuess);
		
		Vector uSmallApproximate = model.solveDirichlet(meshSmall, new Vector2Function(uSmallExtract));
		Tools.plotVector(meshSmall, outputFolder, "u_small_approximate.dat", uSmallApproximate);
		Tools.plotVector(meshSmall, outputFolder, "u_small_guess_approximate_diff.dat", 
				FMath.axpy(-1.0, uSmallApproximate, uSmallGuess));
		
		//TEST 3. Only up side of the domain is Dirichlet boundary
		Function diriBoundaryMark = new AbstractFunction("x","y"){
			@Override
			public double value(Variable v) {
				//double x = v.get("x");
				double y = v.get("y");
				if(Math.abs(y - 3.0) < Constant.eps)
					return 1.0;
				else
					return -1.0;
			}
		};
		Vector uMix = model.solveMixedBorder(meshSmall, 
				diriBoundaryMark, new Vector2Function(uSmallExtract),
				null,FC.c1.D(model.mu_s.M(model.mu_a).M(3.0)));
		Tools.plotVector(meshSmall, outputFolder, "u_mix.dat", uMix);
	}

}
