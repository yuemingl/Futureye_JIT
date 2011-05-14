package edu.uta.futureye.tutorial;

import java.io.File;

import edu.uta.futureye.algebra.SolverJBLAS;
import edu.uta.futureye.algebra.SparseVector;
import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.FX;
import edu.uta.futureye.function.basic.Vector2Function;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.function.intf.VectorFunction;
import edu.uta.futureye.function.operator.FMath;
import edu.uta.futureye.io.MeshReader;
import edu.uta.futureye.io.MeshWriter;
import edu.uta.futureye.lib.assembler.AssemblerScalar;
import edu.uta.futureye.lib.weakform.WeakFormDerivative;
import edu.uta.futureye.util.Utils;
import edu.uta.futureye.util.container.NodeList;

public class Tools {
	public static void plotVector(Mesh mesh, String outputFolder, String fileName, Vector v, Vector ...vs) {
	    MeshWriter writer = new MeshWriter(mesh);
	    if(!outputFolder.isEmpty()) {
		    File file = new File("./"+outputFolder);
			if(!file.exists()) {
				file.mkdirs();
			}
	    }
	    writer.writeTechplot("./"+outputFolder+"/"+fileName, v, vs);
	}

	public static void plotFunction(Mesh mesh, String outputFolder, String fileName, Function fun, Function ...funs) {
	    NodeList list = mesh.getNodeList();
	    int nNode = list.size();
		Variable var = new Variable();
		Vector v = new SparseVector(nNode);
		Vector[] vs = new SparseVector[funs.length];
		for(int i=0;i<funs.length;i++) {
			vs[i] = new SparseVector(nNode);
		}
	    for(int i=1;i<=nNode;i++) {
	    	Node node = list.at(i);
	    	var.setIndex(node.globalIndex);
	    	for(int j=1;j<=node.dim();j++) {
	    		if(fun.varNames().size()==node.dim())
	    			var.set(fun.varNames().get(j-1), node.coord(j));
	    	}
	    	v.set(i, fun.value(var));
	    	for(int j=0;j<funs.length;j++) {
	    		vs[j].set(i,funs[j].value(var));
	    	}
	    }
	    plotVector(mesh,outputFolder,fileName,v,vs);
	}
	
	public static Vector computeDerivative(Mesh mesh, Vector U, String varName) {
		mesh.clearBorderNodeMark();
		
		WeakFormDerivative weakForm = new WeakFormDerivative(varName);
		weakForm.setParam(new Vector2Function(U));
		
		AssemblerScalar assembler = new AssemblerScalar(mesh, weakForm);
		assembler.assemble();
		Matrix stiff = assembler.getStiffnessMatrix();
		Vector load = assembler.getLoadVector();
		
		SolverJBLAS solver = new SolverJBLAS();
		Vector w = solver.solveDGESV(stiff, load);
		return w;
	}
	
	public static Vector computeLaplace2D(Mesh mesh, Vector U) {
		Vector ux = Tools.computeDerivative(mesh,U,"x");
		Vector uy = Tools.computeDerivative(mesh,U,"y");
		Vector uxx = Tools.computeDerivative(mesh,ux,"x");
		Vector uyy = Tools.computeDerivative(mesh,uy,"y");
		Vector LpU = FMath.axpy(1.0, uxx, uyy);
		return LpU;
	}
	
	public static void main(String[] args) {
        MeshReader reader = new MeshReader("triangle.grd");
        Mesh mesh = reader.read2DMesh();
        //fun(x,y)=x^2+y^2
        Function fun = FX.fx.M(FX.fx).A(FX.fy.M(FX.fy));
        plotFunction(mesh,".","testPlotFun.dat",fun);
        VectorFunction gradFun = FMath.grad(fun);
        plotFunction(mesh,".","testPlotFunGrad.dat",
        			gradFun.get(1),gradFun.get(2));
        
        Laplace model = new Laplace();
        model.run();
        Vector u = model.u;
        Vector ux = computeDerivative(model.mesh, u, "x");
        Vector uy = computeDerivative(model.mesh, u, "y");
        Vector uLaplace = computeLaplace2D(model.mesh, u);
        plotVector(model.mesh,".","testPlotUx.dat", ux);
        plotVector(model.mesh,".","testPlotUy.dat", uy);
        plotVector(model.mesh,".","testPlotULaplace.dat", uLaplace);
        //method gaussSmooth() need to know neighbor nodes
        model.mesh.computeNeighborNodes();
        uLaplace = Utils.gaussSmooth(model.mesh, uLaplace, 1, 0.5);
        uLaplace = Utils.gaussSmooth(model.mesh, uLaplace, 1, 0.5);
        uLaplace = Utils.gaussSmooth(model.mesh, uLaplace, 1, 0.5);
        plotVector(model.mesh,".","testPlotULaplaceSmooth.dat", uLaplace);

	}
}
