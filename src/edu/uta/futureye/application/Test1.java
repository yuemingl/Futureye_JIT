package edu.uta.futureye.application;

import java.io.File;
import java.io.FileOutputStream;
import java.io.OutputStreamWriter;
import java.util.HashMap;

import edu.uta.futureye.algebra.Solver;
import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.core.DOF;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.core.NodeType;
import edu.uta.futureye.function.AbstractFunction;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.basic.FDelta;
import edu.uta.futureye.function.basic.FX;
import edu.uta.futureye.function.basic.Vector2Function;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.function.intf.ScalarShapeFunction;
import edu.uta.futureye.function.operator.FOBasic;
import edu.uta.futureye.io.MeshReader;
import edu.uta.futureye.io.MeshWriter;
import edu.uta.futureye.lib.assembler.AssemblerScalar;
import edu.uta.futureye.lib.shapefun.SFLinearLocal2D;
import edu.uta.futureye.lib.weakform.WeakFormL22D;
import edu.uta.futureye.lib.weakform.WeakFormLaplace2D;
import edu.uta.futureye.util.container.NodeList;

public class Test1 {
	public static void allTest() {
		MeshReader reader = new MeshReader("prostate_test1.grd");
		Mesh mesh = reader.read2DMesh();
		mesh.computeNodeBelongsToElements();
		
		//Mark border type
		HashMap<NodeType, Function> mapNTF = new HashMap<NodeType, Function>();
		mapNTF.put(NodeType.Robin, null);
		
		//mapNTF.put(NodeType.Dirichlet, null);		
		
//		mapNTF.put(NodeType.Dirichlet, new FAbstract("x","y"){
//			@Override
//			public double value(Variable v) {
//				//if(Math.abs(v.get("y"))<0.01 || Math.abs(v.get("x")-5.0)<0.01)
//					if(Math.abs(v.get("y"))<0.01)
//					return 1.0;
//				else
//					return -1.0;
//			}
//		});
		
		
		mesh.markBorderNode(mapNTF);

		SFLinearLocal2D[] shapeFun = new SFLinearLocal2D[3];
		shapeFun[0] = new SFLinearLocal2D(1);
		shapeFun[1] = new SFLinearLocal2D(2);
		shapeFun[2] = new SFLinearLocal2D(3);
		
		//Asign degree of freedom to element
		for(int i=1;i<=mesh.getElementList().size();i++) {
			Element e = mesh.getElementList().at(i);
			for(int j=1;j<=e.nodes.size();j++) {
				//Asign shape function to DOF
				DOF dof = new DOF(j,e.nodes.at(j).globalIndex,shapeFun[j-1]);
				e.addNodeDOF(j, dof);
			}
		}
		
		//User defined weak form of PDE (including bounder conditions)
		WeakFormLaplace2D weakForm = new WeakFormLaplace2D();
		//Right hand side
		Variable x0 = new Variable();
		x0.set("x", 1.0);
		x0.set("y", 2.8);
		FDelta delta = new FDelta(x0,0.01,2e5);
		weakForm.setF(delta);
		
		final double cx[] = {2.0, 3.0, 4.0};
		double cy[] = {1.5, 2.0, 2.5};
		double mu_a[] = {0.1,0.2,0.3,0.4,1.0};
		for(int cxi=0;cxi<cx.length;cxi++)
			for(int cyi=0;cyi<cy.length;cyi++)
				for(int mu_ai=0;mu_ai<mu_a.length;mu_ai++) {
			
			final double fcx = cx[cxi];
			final double fcy = cy[cyi];
			final double fmu_a = mu_a[mu_ai];
									
			weakForm.setParam(
			new FC(0.02), //  1/(3*mu_s') = 0.02
			new AbstractFunction("x","y"){ //mu_a
				@Override
				public double value(Variable v) {
					double dx = v.get("x")-fcx;
					double dy = v.get("y")-fcy;
					if(Math.sqrt(dx*dx+dy*dy)<0.5)
						return fmu_a; 
					else
						return 0.1;
				}
			},
			new FC(0.05),null // d*u + k*u_n= q (自然边界：d==k, q=0)
			); 
			
			AssemblerScalar assembler = new AssemblerScalar(mesh, weakForm);
			System.out.println("Begin Assemble...");
			Matrix stiff = assembler.getStiffnessMatrix();
			Vector load = assembler.getLoadVector();
			assembler.imposeDirichletCondition(new FC(0.0));
			System.out.println("Assemble done!");
			
			Solver solver = new Solver();
	//		Matrix stiff2 = new Matrix(stiff.getRowDim()+1,stiff.getColDim()+1);
	//		Vector load2 = new Vector(load.getDim()+1);
	//		for(int i=1;i<=stiff.getColDim();i++) {
	//			for(int j=1;j<=stiff.getRowDim();j++)
	//				stiff2.set(i, j, stiff.get(i, j));
	//			load2.set(i, load.get(i));
	//		}
	//		for(int i=1;i<=stiff2.getColDim();i++) {
	//			stiff2.set(stiff2.getRowDim(), i, 1.0);
	//			stiff2.set(i, stiff2.getColDim(), 1.0);
	//		}
	//		load2.set(load2.getDim(), 1000000);
	//		Vector u = solver.solve(stiff2, load2);
			
			
	//		for(int i=1;i<=stiff.getColDim();i++) {
	//			double a = stiff.get(1, i);
	//			stiff.plusValue(2, i, a);
	//			stiff.set(1, i, 1.0);
	//		}
	//		load.plusValue(2, load.get(1));
	//		load.set(1, 1000);
			
			Vector u = solver.solve(stiff, load);
		    
			System.out.println("u=");
		    for(int i=1;i<=u.getDim();i++)
		        System.out.println(String.format("%.3f", u.get(i)));	
		    
		    MeshWriter writer = new MeshWriter(mesh);
		    writer.writeTechplot("prostate_test1"+String.format("_%d_%d_%d", cxi,cyi,mu_ai)+".dat", 
		    		u);
		}
		
	}

	
	public static void paramInverseTest() {
		MeshReader reader = new MeshReader("prostate_test1.grd");
		Mesh mesh = reader.read2DMesh();
		mesh.computeNodeBelongsToElements();
		
		//Mark border type
		HashMap<NodeType, Function> mapNTF = new HashMap<NodeType, Function>();
		mapNTF.put(NodeType.Robin, null);
		mesh.markBorderNode(mapNTF);

		SFLinearLocal2D[] shapeFun = new SFLinearLocal2D[3];
		shapeFun[0] = new SFLinearLocal2D(1);
		shapeFun[1] = new SFLinearLocal2D(2);
		shapeFun[2] = new SFLinearLocal2D(3);
		
		//Asign degree of freedom to element
		for(int i=1;i<=mesh.getElementList().size();i++) {
			Element e = mesh.getElementList().at(i);
			for(int j=1;j<=e.nodes.size();j++) {
				//Asign shape function to DOF
				DOF dof = new DOF(j,e.nodes.at(j).globalIndex,shapeFun[j-1]);
				e.addNodeDOF(j, dof);
			}
		}
		
		//User defined weak form of PDE (including bounder conditions)
		WeakFormLaplace2D weakForm = new WeakFormLaplace2D();
		//Right hand side
		Variable x0 = new Variable();
		x0.set("x", 1.0);
		x0.set("y", 2.8);
		FDelta delta = new FDelta(x0,0.01,2e5);
		weakForm.setF(delta);
		
		final double fcx = 2.0;
		final double fcy = 2.5;
		final double fmu_a = 1.0;

		weakForm.setParam(
		new FC(0.02), //  1/(3*mu_s') = 0.02
		new AbstractFunction("x","y"){ //mu_a
			@Override
			public double value(Variable v) {
				double dx = v.get("x")-fcx;
				double dy = v.get("y")-fcy;
				if(Math.sqrt(dx*dx+dy*dy)<0.5)
					return fmu_a; 
				else
					return 0.1;
			}
		},
		new FC(0.05),null // d*u + k*u_n= q (自然边界：d==k, q=0)
		); 
		
		AssemblerScalar assembler = new AssemblerScalar(mesh, weakForm);
		System.out.println("Begin Assemble...");
		Matrix stiff = assembler.getStiffnessMatrix();
		Vector load = assembler.getLoadVector();
		assembler.imposeDirichletCondition(new FC(0.0));
		System.out.println("Assemble done!");
		
		Solver solver = new Solver();
		Vector u = solver.solve(stiff, load);
		System.out.println("u=");
	    for(int i=1;i<=u.getDim();i++)
	        System.out.println(String.format("%.3f", u.get(i)));	
	    
	    MeshWriter writer = new MeshWriter(mesh);
	    writer.writeTechplot("prostate_test1"+ String.format("_u.dat"), u);
	    
	    
		//User defined weak form of PDE (including bounder conditions)
		HashMap<NodeType, Function> mapNTF2 = new HashMap<NodeType, Function>();
		mapNTF2.put(NodeType.Dirichlet, null);
		mesh.clearBorderNodeMark();
		mesh.markBorderNode(mapNTF2);
		
		WeakFormL22D weakFormL2 = new WeakFormL22D();
		weakFormL2.setF(delta);
		weakFormL2.setParam(
				new FC(0.02), 
				new Vector2Function(u)
				);
		
		AssemblerScalar assembler2 = new AssemblerScalar(mesh, weakFormL2);
		System.out.println("Begin Assemble...");
		Matrix stiff2 = assembler2.getStiffnessMatrix();
		Vector load2 = assembler2.getLoadVector();
		assembler2.imposeDirichletCondition(new FC(0.1));
		System.out.println("Assemble done!");
		
		Solver solver2 = new Solver();
		Vector u2 = solver2.solve(stiff2, load2);
		System.out.println("alpha=");
	    for(int i=1;i<=u2.getDim();i++)
	        System.out.println(String.format("%.3f", u2.get(i)));	

	    writer.writeTechplot("prostate_test1"+String.format("_alpha.dat"), u2);
	}
	
	public static void testWeakFormGCM() {
		MeshReader reader = new MeshReader("triangle.grd");
		Mesh mesh = reader.read2DMesh();
		mesh.computeNodeBelongsToElements();
		
		HashMap<NodeType, Function> mapNTF = new HashMap<NodeType, Function>();
		mapNTF.put(NodeType.Dirichlet, null);		
		mesh.markBorderNode(mapNTF);

		SFLinearLocal2D[] shapeFun = new SFLinearLocal2D[3];
		for(int i=0;i<3;i++)
			shapeFun[i] = new SFLinearLocal2D(i+1);
		
		//Asign degree of freedom to element
		for(int i=1;i<=mesh.getElementList().size();i++) {
			Element e = mesh.getElementList().at(i);
			for(int j=1;j<=e.nodes.size();j++) {
				//Asign shape function to DOF
				DOF dof = new DOF(j,e.nodes.at(j).globalIndex,shapeFun[j-1]);
				e.addNodeDOF(j, dof);
			}
		}
		
		Function x = new FX("x");
		Function y = new FX("y");
		Function x2 = FOBasic.Mult(x,x);
		Function y2 = FOBasic.Mult(y,y);
		
		Function f =
				FOBasic.PlusAll(
						FOBasic.Mult(new FC(2.0), x2),
						FOBasic.Mult(new FC(2.0), y2),
						FOBasic.Mult(new FC(2.0), FOBasic.Mult(x2,y)),
						FOBasic.Mult(new FC(2.0), FOBasic.Mult(y2,x)),
						FOBasic.Mult(new FC(-18.0), x),
						FOBasic.Mult(new FC(-18.0), y),
						new FC(-36.0)
				);
		System.out.println(f);
		
		//User defined weak form of PDE (including bounder conditions)
		//-k*\Delta{u} + b\dot\Nabla{u} + c*u = f
		//u(x,y)=0, (x,y)\in\partial{\Omega}
		//k=-1
		//b=(1 1)'
		//c=0
		//f=2*(x^2+y^2) + (2xy-18)*(x+y) -36
		//u=(x^2-9)*(y^2-9)
		WeakFormGCM weakForm = new WeakFormGCM();
		weakForm.setF(f);
		weakForm.setParam(
				new FC(-1.0),//注意，有负号!!!
				new FC(0.0),
				new FC(1.0),
				new FC(1.0)
			);
		
		mesh.clearBorderNodeMark();
		mesh.markBorderNode(mapNTF);

		AssemblerScalar assembler = new AssemblerScalar(mesh, weakForm);
		System.out.println("Begin Assemble...solveGCM");
		Matrix stiff = assembler.getStiffnessMatrix();
		Vector load = assembler.getLoadVector();
		//Dirichlet condition
		assembler.imposeDirichletCondition(new FC(0.0));
		System.out.println("Assemble done!");

		Solver solver = new Solver();
		Vector u = solver.solve(stiff, load);
		
	    MeshWriter writer = new MeshWriter(mesh);
	    writer.writeTechplot("testWeakFormGCM.dat", u);
	}
	
	
	public static void testLeastSquare() {
		String gridFile = "prostate_test_least_square.grd";
		
		MeshReader reader = null;
		reader = new MeshReader(gridFile);
		Mesh mesh = reader.read2DMesh();
		
		ScalarShapeFunction[] shapeFun = null;
		mesh.computeNodeBelongsToElements();
		mesh.computeNeighborNodes();
		//Assign degree of freedom to element
		shapeFun = new SFLinearLocal2D[3];
		for(int i=0;i<3;i++)
			shapeFun[i] = new SFLinearLocal2D(i+1);
		
		//Assign shape function to DOF
		for(int i=1;i<=mesh.getElementList().size();i++) {
			Element e = mesh.getElementList().at(i);
			for(int j=1;j<=e.nodes.size();j++) {
				DOF dof = new DOF(j,e.nodes.at(j).globalIndex,shapeFun[j-1]);
				e.addNodeDOF(j, dof);
			}
		}

		Model model = new Model();
		model.outputFolder = "prostate_least_square";
	
		NodeList list =mesh.getNodeList();
		int nNode = list.size();
		StringBuilder nd = new StringBuilder();
		for(int i=1;i<=nNode;i++) {
			Node node = list.at(i);
			nd.append(String.format("%d\t%f\t%f\r\n", node.globalIndex,
						node.coord(1),node.coord(2)));
		}
		writeStringToFile(model.outputFolder+"\\node.txt",nd.toString());
		
		StringBuilder rhs = new StringBuilder();
		StringBuilder bkU = new StringBuilder();
		rhs.append("  ");
		bkU.append("  ");
		
		double deta = 0.1;
		for(int l=1;l<=15;l++) {
			
			//----------------Left light source ---------------------
			model.setDelta(1.0+deta*l, 2.8);
			//Solve background forward problem
			model.setMu_a(0.0, 0.0, 0.0, 
					0.1, 1);
			Vector bkUL = model.solveForwardNeumann(mesh);
			model.plotVector(mesh, bkUL, "bkUL"+l+".dat");
			
			//Solve forward problem with inclusion
			model.setMu_a(2.6, 2.7, 0.2, 
					2.0, 1);
			if(l==1)
				model.plotFunction(mesh, model.mu_a, "alpha_real.dat");
			Vector incUL = model.solveForwardNeumann(mesh);
			model.plotVector(mesh, incUL, "incUL"+l+".dat");
			
			rhs.delete(0, rhs.length()-1);
			bkU.delete(0, bkU.length()-1);
			for(int i=1;i<=nNode;i++) {
				Node node = list.at(i);
				if(Math.abs(3.0-node.coord(2))<0.01)
					rhs.append(String.format("%d\t%f\r\n", node.globalIndex,
							incUL.get(i)-bkUL.get(i)));
				bkU.append(String.format("%d\t%f\r\n",node.globalIndex, bkUL.get(i)));
			}
			writeStringToFile(model.outputFolder+"\\rhsL"+l+".txt",rhs.toString());
			writeStringToFile(model.outputFolder+"\\bkUL"+l+".txt",bkU.toString());

		
			//----------------Right light source ---------------------
			model.setDelta(4.0-deta*l, 2.8);
			
			model.setMu_a(0.0, 0.0, 0.0, 
					0.1, 1);
			Vector bkUR = model.solveForwardNeumann(mesh);
			model.plotVector(mesh, bkUR, "bkUR"+l+".dat");
			
			model.setMu_a(2.6, 2.7, 0.2, 
					2.0, 1);
			Vector incUR = model.solveForwardNeumann(mesh);
			model.plotVector(mesh, incUR, "incUR"+l+".dat");
			
			rhs.delete(0, rhs.length()-1);
			bkU.delete(0, bkU.length()-1);
			for(int i=1;i<=nNode;i++) {
				Node node = list.at(i);
				if(Math.abs(3.0-node.coord(2))<0.01)
					rhs.append(String.format("%d\t%f\r\n", node.globalIndex,
							incUR.get(i)-bkUR.get(i)));
				bkU.append(String.format("%d\t%f\r\n",node.globalIndex, bkUR.get(i)));
			}
			writeStringToFile(model.outputFolder+"\\rhsR"+l+".txt",rhs.toString());
			writeStringToFile(model.outputFolder+"\\bkUR"+l+".txt",bkU.toString());
		}
	}
	
	/**
	 * one direction
	 */
	public static void testLeastSquare2() {
		String gridFile = "prostate_test_least_square.grd";
		//String gridFile = "prostate_test3.grd";
		
		MeshReader reader = null;
		reader = new MeshReader(gridFile);
		Mesh mesh = reader.read2DMesh();
		
		ScalarShapeFunction[] shapeFun = null;
		mesh.computeNodeBelongsToElements();
		mesh.computeNeighborNodes();
		//Assign degree of freedom to element
		shapeFun = new SFLinearLocal2D[3];
		for(int i=0;i<3;i++)
			shapeFun[i] = new SFLinearLocal2D(i+1);
		
		//Assign shape function to DOF
		for(int i=1;i<=mesh.getElementList().size();i++) {
			Element e = mesh.getElementList().at(i);
			for(int j=1;j<=e.nodes.size();j++) {
				DOF dof = new DOF(j,e.nodes.at(j).globalIndex,shapeFun[j-1]);
				e.addNodeDOF(j, dof);
			}
		}

		Model model = new Model();
		model.outputFolder = "prostate_least_square";
	
		NodeList list =mesh.getNodeList();
		int nNode = list.size();
		StringBuilder nd = new StringBuilder();
		for(int i=1;i<=nNode;i++) {
			Node node = list.at(i);
			nd.append(String.format("%d\t%f\t%f\r\n", node.globalIndex,
						node.coord(1),node.coord(2)));
		}
		writeStringToFile(model.outputFolder+"\\node.txt",nd.toString());
		
		StringBuilder rhs = new StringBuilder();
		StringBuilder bkU = new StringBuilder();
		rhs.append("  ");
		bkU.append("  ");
		
		double deta = 0.2;
		for(int l=1;l<=5;l++) {
			
			//----------------Left light source ---------------------
			model.setDelta(1.0+deta*l, 2.8);
			//Solve background forward problem
			model.setMu_a(0.0, 0.0, 0.0, 
					0.1, 1);
			Vector bkUL = model.solveForwardNeumann(mesh);
			model.plotVector(mesh, bkUL, "bkUL"+l+".dat");
			
			//Solve forward problem with inclusion
			model.setMu_a(3.0, 2.7, 0.2, 
					2.0, 1);
			if(l==1)
				model.plotFunction(mesh, model.mu_a, "alpha_real.dat");
			Vector incUL = model.solveForwardNeumann(mesh);
			model.plotVector(mesh, incUL, "incUL"+l+".dat");
			
			rhs.delete(0, rhs.length()-1);
			bkU.delete(0, bkU.length()-1);
			for(int i=1;i<=nNode;i++) {
				Node node = list.at(i);
				if(Math.abs(3.0-node.coord(2))<0.01)
					rhs.append(String.format("%d\t%f\r\n", node.globalIndex,
							incUL.get(i)-bkUL.get(i)));
				bkU.append(String.format("%d\t%f\r\n",node.globalIndex, bkUL.get(i)));
			}
			writeStringToFile(model.outputFolder+"\\rhsL"+l+".txt",rhs.toString());
			writeStringToFile(model.outputFolder+"\\bkUL"+l+".txt",bkU.toString());
		
		}
	}

	public static void writeStringToFile(String fileName,String content) {
		FileOutputStream out;
		try {
			out = new FileOutputStream(new File(fileName));
			OutputStreamWriter writer = new OutputStreamWriter(out, "UTF-8");
			writer.write(content);
			writer.close();
			out.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public static void main(String[] args) {
		
		//testWeakFormGCM();
		
		Model model = new Model();
				
		//3*5 => 30*50
		//model.run(1,1,"prostate_test1");
		//3*5 => 10*15
		//model.run(1,2,"prostate_test2");
		//3*5 => 15*25
		//model.run(1,3,"prostate_test3");
		//3*5 => manual adaptive
		//model.run(1,4,"prostate_test4_linear");
		//model.run(2,4,"prostate_test4_quadratic");
		//3*5 => mixed
		//model.run(1,5,"prostate_test5_mixed");
		//3*5 => rectangle
		//model.run(1,6,"prostate_test6_rectangle");
		
		
		//model.runAdaptive(1,1,"prostate_test1");
		//model.runAdaptive(1,2,"prostate_test2");
		//model.runAdaptive(1,3,"prostate_test3");
		//model.runAdaptive(1,2,"prostate_test2");
		//model.runAdaptive(1,6,"prostate_test6_rectangle");
		//model.runAdaptive(1,7,"prostate_test7_rectangle");
		
		//model.runAdaptive(1,5,"prostate_test5_mixed");
	
		model.runGCMTest(1,
				"prostate_test3_ex.grd",
				"prostate_test3.grd",
				"prostate_test3_gcm");
		
		
		//testLeastSquare();
		//testLeastSquare2();
	}
	
	
	
}
