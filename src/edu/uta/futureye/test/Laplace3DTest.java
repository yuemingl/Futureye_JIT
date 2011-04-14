package edu.uta.futureye.test;

import java.util.HashMap;

import edu.uta.futureye.algebra.CompressedRowMatrix;
import edu.uta.futureye.algebra.FullVector;
import edu.uta.futureye.algebra.Solver;
import edu.uta.futureye.algebra.SparseMatrix;
import edu.uta.futureye.algebra.SparseVector;
import edu.uta.futureye.algebra.intf.AlgebraMatrix;
import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.core.DOF;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.NodeType;
import edu.uta.futureye.core.intf.Assembler;
import edu.uta.futureye.function.AbstractFunction;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.io.MeshReader;
import edu.uta.futureye.io.MeshWriter;
import edu.uta.futureye.lib.assembler.AssemblerScalar;
import edu.uta.futureye.lib.element.FELinearTetrahedron;
import edu.uta.futureye.lib.shapefun.SFLinearLocal3D;
import edu.uta.futureye.lib.weakform.WeakFormLaplace;
import edu.uta.futureye.lib.weakform.WeakFormLaplace3D;
import edu.uta.futureye.util.container.ElementList;
import edu.uta.futureye.util.container.NodeList;

public class Laplace3DTest {
	
	public static void triangleTest() {
		String meshName = "block1";
		MeshReader reader = new MeshReader(meshName+".grd");
		Mesh mesh = reader.read3DMesh(); //3D
		mesh.computeNodeBelongsToElements(); //worked in 3D
		
		HashMap<NodeType, Function> mapNTF = new HashMap<NodeType, Function>();
		mapNTF.put(NodeType.Robin, new AbstractFunction("x","y","z"){
			@Override
			public double value(Variable v) {
				if(1.0-v.get("x")<0.01)
					return 1.0;
				else
					return -1.0;
			}
		});
		mapNTF.put(NodeType.Dirichlet, null);
		mesh.markBorderNode(mapNTF);
		
		NodeList nList = mesh.getNodeList();
		for(int i=1;i<=nList.size();i++) {
			if(nList.at(i).getNodeType() == NodeType.Robin)
				System.out.println(nList.at(i));
		}

//		//Asign degree of freedom to element
//		SFLinearLocal3D[] shapeFun = new SFLinearLocal3D[4];
//		for(int i=0;i<4;i++)
//			shapeFun[i] = new SFLinearLocal3D(i+1);
//		for(int i=1;i<=mesh.getElementList().size();i++) {
//			Element e = mesh.getElementList().at(i);
//			for(int j=1;j<=e.nodes.size();j++) {
//				//Asign shape function to DOF
//				DOF dof = new DOF(j,e.nodes.at(j).globalIndex,shapeFun[j-1]);
//				e.addNodeDOF(j, dof);
//			}
//		}
		
		//Use element library
		ElementList eList = mesh.getElementList();
		FELinearTetrahedron fe = new FELinearTetrahedron();
		for(int i=1;i<=eList.size();i++)
			fe.assignTo(eList.at(i));
		
		//User defined weak form of PDE (including bounder conditions)
		//-\Delta{u} = f
		//u(x,y)=0, (x,y)\in\partial{\Omega}
		//u=(x^2-9)*(y^2-9)
		//f=-2*(x^2+y^2)+36
		WeakFormLaplace3D weakForm = new WeakFormLaplace3D();
		
		weakForm.setF(new FC(1.0));
		
		weakForm.setParam(
					new FC(1.0),
					new FC(1.0),
					new FC(0.0),
					null //Robin: 6*y^2-54
				);
		
		Assembler assembler = new AssemblerScalar(mesh, weakForm);
		System.out.println("Assemble...");
		long begin = System.currentTimeMillis();
		assembler.assemble();
		Matrix stiff = assembler.getStiffnessMatrix();
		Vector load = assembler.getLoadVector();
		long end = System.currentTimeMillis();
		System.out.println("Assemble done!");
		System.out.println("Assemble time:"+(end-begin));
		
		System.out.println("Impose Dirichlet condition...");
		assembler.imposeDirichletCondition(new FC(0.0));
		System.out.println("Impose Dirichlet condition done!");

		//Initial value for iteration solvers
		SparseVector u = new SparseVector(load.getDim(),0.0003);
		
		
		Solver solver = new Solver();
		begin = System.currentTimeMillis();
		
		//CG
		System.out.println("begin construct AlgebraMatrix...");
		AlgebraMatrix algStiff = new CompressedRowMatrix((SparseMatrix)stiff,true);
		System.out.println("end construct AlgebraMatrix!");
		FullVector algLoad = new FullVector((SparseVector)load);
		FullVector algU = new FullVector(u);
		solver.CG(algStiff, algLoad, algU);
		double[] data = algU.getData();
		for(int i=0;i<data.length;i++) {
			u.set(i+1, data[i]);
		}
		
		//CG1
		//solver.CG1(stiff, load, u);
		
		//CG2
		//solver.CG2((SparseMatrix)stiff, (SparseVector)load, (SparseVector)u);
		
		//Blas Full
		//u = (SparseVector) solver.solve(stiff, load);
		
		
		end = System.currentTimeMillis();
		//Solve time used:5214
		System.out.println("Solve time used:"+(end-begin));
		

	    System.out.println("u=");
	    for(int i=1;i<=u.getDim();i++)
	        System.out.println(String.format("%.3f", u.get(i)));	
	    
	    MeshWriter writer = new MeshWriter(mesh);
	    writer.writeTechplot(meshName+"_out.dat", u);
		
	}
	
	public static void triangleTest_WeakFormxD() {
		String meshName = "block1";
		MeshReader reader = new MeshReader(meshName+".grd");
		Mesh mesh = reader.read3DMesh(); //3D
		mesh.computeNodeBelongsToElements(); //worked in 3D
		
		HashMap<NodeType, Function> mapNTF = new HashMap<NodeType, Function>();
		mapNTF.put(NodeType.Dirichlet, null);
		mesh.markBorderNode(mapNTF);
		
//		NodeList nList = mesh.getNodeList();
//		for(int i=1;i<=nList.size();i++) {
//			if(nList.at(i).getNodeType() == NodeType.Dirichlet)
//				System.out.println(nList.at(i));
//		}

		SFLinearLocal3D[] shapeFun = new SFLinearLocal3D[4];
		for(int i=0;i<4;i++)
			shapeFun[i] = new SFLinearLocal3D(i+1);
		
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
		//-\Delta{u} = f
		//u(x,y)=0, (x,y)\in\partial{\Omega}
		//u=(x^2-9)*(y^2-9)
		//f=-2*(x^2+y^2)+36
		WeakFormLaplace weakForm = new WeakFormLaplace();
		
		weakForm.setF(new FC(1.0));
		
		weakForm.setParam(
					null,
					new FC(1.0),
					new FC(0.0),
					null //Robin: 6*y^2-54
				);
		
		
		AssemblerScalar assembler = new AssemblerScalar(mesh, weakForm);
		System.out.println("Begin Assemble...");
		long begin = System.currentTimeMillis();
		assembler.assemble();
		Matrix stiff = assembler.getStiffnessMatrix();
		Vector load = assembler.getLoadVector();
		assembler.imposeDirichletCondition(new FC(0.0));
		long end = System.currentTimeMillis();
		System.out.println("Assemble done!");
		System.out.println("Time used:"+(end-begin));
		
		Solver solver = new Solver();
		Vector u = solver.solve(stiff, load);
	    System.out.println("u=");
	    for(int i=1;i<=u.getDim();i++)
	        System.out.println(String.format("%.3f", u.get(i)));	
	    
	    MeshWriter writer = new MeshWriter(mesh);
	    writer.writeTechplot(meshName+"_out.dat", u);
		
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		//triangleTest();
		triangleTest_WeakFormxD();
	}

}
