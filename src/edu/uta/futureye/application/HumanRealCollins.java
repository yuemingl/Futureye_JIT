package edu.uta.futureye.application;

import static edu.uta.futureye.function.FMath.C;
import static edu.uta.futureye.function.FMath.C0;

import java.util.HashMap;

import edu.uta.futureye.algebra.SparseVectorHashMap;
import edu.uta.futureye.algebra.intf.SparseMatrix;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.algebra.solver.Solver;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.NodeType;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.basic.FDelta;
import edu.uta.futureye.function.basic.Vector2MathFunc;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.io.MeshReader;
import edu.uta.futureye.io.MeshWriter;
import edu.uta.futureye.lib.assembler.AssemblerScalar;
import edu.uta.futureye.lib.element.FEBilinearRectangle;
import edu.uta.futureye.lib.element.FELinearTetrahedron;
import edu.uta.futureye.lib.weakform.WeakFormLaplace2D;
import edu.uta.futureye.lib.weakform.WeakFormLaplace3D;
import edu.uta.futureye.util.container.ElementList;

public class HumanRealCollins {
	public static String folder = "HumanRealCollins";
	public static String outputFolder = folder+"/output";

	public Mesh read3DMesh(String file) {
		MeshReader reader = new MeshReader(folder+"/"+file);
		Mesh mesh = reader.read3DMesh();
		Vector v = new SparseVectorHashMap(mesh.getNodeList().size());
		Tools.plotVector(mesh, outputFolder, "iso2mesh_head.dat", v);
		return mesh;
	}
	
	
	public void test1(Mesh mesh) {
        //Compute geometry relationship of nodes and elements
        mesh.computeNodeBelongsToElements();

        //Mark border types
        HashMap<NodeType, MathFunc> mapNTF =
                new HashMap<NodeType, MathFunc>();
        mapNTF.put(NodeType.Dirichlet, null);
        mesh.markBorderNode(mapNTF);

        //Use element library to assign degrees of
        //  freedom (DOF) to element
		ElementList eList = mesh.getElementList();
		FELinearTetrahedron fe = new FELinearTetrahedron();
		for(int i=1;i<=eList.size();i++)
			fe.assignTo(eList.at(i));
		
		
		//User defined weak form of PDE (including bounder conditions)
		WeakFormLaplace3D weakForm = new WeakFormLaplace3D();
		
		//Function f = new FDelta(new Variable("x",90.9314).set("y",230.306).set("z",59.5004),1e-2,1e5);
		MathFunc f = new FDelta(new Variable("x",105.572).set("y",52.3113).set("z",152.091),1e1,1e5);
		Tools.plotFunction(mesh, outputFolder, "lightsource.dat", f);
		weakForm.setF(f);
		
		weakForm.setParam(
					new FC(20),
					new FC(0.1),
					new FC(0.0),
					new FC(20)
				);

        //Assembly process
        AssemblerScalar assembler =
                new AssemblerScalar(mesh, weakForm);
        assembler.assemble();
        SparseMatrix stiff = assembler.getStiffnessMatrix();
        Vector load = assembler.getLoadVector();
        //Boundary condition
        assembler.imposeDirichletCondition(C0);

        //Solve linear system
        Solver solver = new Solver();
        Vector u = solver.solveCG(stiff, load);
        System.out.println("u=");
        for(int i=1;i<=u.getDim();i++)
            System.out.println(String.format("%.3f", u.get(i)));

        //Output results to an Techplot format file
		Tools.plotVector(mesh, outputFolder, "forward.dat", u);

	}
	
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		HumanRealCollins a = new HumanRealCollins();
		//a.read3DMesh("collins_brain.grd");
		Mesh mesh = a.read3DMesh("iso2mesh_headv2.grd");
		a.test1(mesh);

	}

}
