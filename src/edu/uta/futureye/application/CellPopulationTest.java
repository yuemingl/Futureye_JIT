package edu.uta.futureye.application;

import static edu.uta.futureye.function.operator.FMath.C0;
import static edu.uta.futureye.function.operator.FMath.X;
import static edu.uta.futureye.function.operator.FMath.Y;
import static edu.uta.futureye.function.operator.FMath.C;

import java.util.HashMap;


import edu.uta.futureye.algebra.SparseVectorHashMap;
import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.SparseMatrix;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.algebra.solver.Solver;
import edu.uta.futureye.algebra.solver.external.SolverJBLAS;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.core.NodeType;
import edu.uta.futureye.function.AbstractMathFunc;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.Vector2Function;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.io.MatlabMatFileReader;
import edu.uta.futureye.io.MeshReader;
import edu.uta.futureye.io.MeshWriter;
import edu.uta.futureye.lib.assembler.AssemblerScalar;
import edu.uta.futureye.lib.element.FEBilinearRectangle;
import edu.uta.futureye.lib.element.FELinearTriangle;
import edu.uta.futureye.lib.weakform.WeakFormLaplace2D;
import edu.uta.futureye.tutorial.T02Laplace;
import edu.uta.futureye.util.Constant;
import edu.uta.futureye.util.container.ElementList;

public class CellPopulationTest {
	
	/**
	 * Read in the matrix of lambda (a picture of cells population)
	 * Elements: 70*287
	 * Nodes: 71*288
	 * @param mesh
	 * @return
	 */
	public static Vector readLambda(Mesh mesh) {
    	MatlabMatFileReader mr = new MatlabMatFileReader("./CellPopulation/lambda-Intensity-71by288.mat");
    	Matrix lambda = mr.getMatrix("I");
    	int m = lambda.getRowDim();
    	int n = lambda.getColDim();
    	Vector v = new SparseVectorHashMap(m*n);
    	double dx = 288.0/287.0;
    	double dy = 71.0/70.0;
    	for(int i=1;i<=m;i++) {
    		for(int j=1;j<=n;j++) {
    			double[] coord = {(j-1)*dx, (i-1)*dy};
    			Node node = mesh.findNodeByCoord(coord, 0.5);
    			if(node == null) {
    				System.out.println(coord[0]+" "+coord[1]);
    				return null;
    			}
    			v.set(node.getIndex(),lambda.get(m-i+1, j)+128.0); //上下颠倒
    		}
    	}
    	return v;
	}
		
	public static void test1() {
        //1.Read in a triangle mesh from an input file with
        //  format ASCII UCD generated by Gridgen
        MeshReader reader = new MeshReader("./CellPopulation/test1.grd");
        Mesh mesh = reader.read2DMesh();
        //Compute geometry relationship of nodes and elements
        mesh.computeNodeBelongsToElements();

        //2.Mark border types
        HashMap<NodeType, MathFunc> mapNTF =
                new HashMap<NodeType, MathFunc>();
//        mapNTF.put(NodeType.Dirichlet, new AbstractFunction("x","y") {
//        	public double value(Variable v) {
//        		double x = v.get("x");
//        		//double y = v.get("y");
//        		if(Math.abs(x-288)<Constant.meshEps)
//        			return 1;
//        		else return 0;
//        	}
//        });
        mapNTF.put(NodeType.Neumann, null);
        mesh.markBorderNode(mapNTF);

        //3.Use element library to assign degrees of
        //  freedom (DOF) to element
        ElementList eList = mesh.getElementList();
        FEBilinearRectangle feLT = new FEBilinearRectangle();
        for(int i=1;i<=eList.size();i++)
            feLT.assignTo(eList.at(i));

        //4.Weak form
        WeakFormLaplace2D weakForm = new WeakFormLaplace2D();
        //参数很重要，否则计算结果不正确
        weakForm.setParam(C(0.02), C(0.1), null, null);
        //Right hand side(RHS): f = lambda
        Vector vf = readLambda(mesh);
        Tools.plotVector(mesh, "./CellPopulation", "f.dat", vf);
        weakForm.setF(new Vector2Function(vf));

        //5.Assembly process
        AssemblerScalar assembler =
                new AssemblerScalar(mesh, weakForm);
        assembler.assemble();
        SparseMatrix stiff = assembler.getStiffnessMatrix();
        Vector load = assembler.getLoadVector();
        //Boundary condition
        assembler.imposeDirichletCondition(C0);

        //6.Solve linear system
        Solver solver = new Solver();
        Vector u = solver.solveCG(stiff, load);
        System.out.println("u=");
        for(int i=1;i<=u.getDim();i++)
            System.out.println(String.format("%.3f", u.get(i)));

        //7.Output results to an Techplot format file
        MeshWriter writer = new MeshWriter(mesh);
        writer.writeTechplot("./CellPopulation/test1.dat", u);

	}
	
    public static void main(String[] args) {
    	test1();

    }
}
