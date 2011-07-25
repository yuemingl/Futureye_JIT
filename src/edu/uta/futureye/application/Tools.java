package edu.uta.futureye.application;

import java.io.File;

import edu.uta.futureye.algebra.Solver;
import edu.uta.futureye.algebra.SolverJBLAS;
import edu.uta.futureye.algebra.SparseVector;
import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.Vector2Function;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.function.operator.FMath;
import edu.uta.futureye.io.MeshWriter;
import edu.uta.futureye.lib.assembler.AssemblerScalar;
import edu.uta.futureye.lib.weakform.WeakFormDerivative;
import edu.uta.futureye.util.container.NodeList;

public class Tools {
	public static Vector computeDerivative(Mesh mesh, Vector U, String varName) {
		//2011-06-29
		//没必要清楚边界标记
		//mesh.clearBorderNodeMark();
		
		WeakFormDerivative weakForm = new WeakFormDerivative(varName);
		weakForm.setParam(new Vector2Function(U));
		
		AssemblerScalar assembler = new AssemblerScalar(mesh, weakForm);
		assembler.assemble();
		Matrix stiff = assembler.getStiffnessMatrix();
		Vector load = assembler.getLoadVector();
		
		//Solver solver = new Solver();
		//Vector w = solver.solveCGS(stiff, load);
		
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
	    for(int i=1;i<=nNode;i++) {
	    	Node node = list.at(i);
	    	var.setIndex(node.globalIndex);
	    	for(int j=1;j<=node.dim();j++) {
	    		if(fun.varNames().size()==node.dim())
	    			var.set(fun.varNames().get(j-1), node.coord(j));
	    	}
	    	v.set(i, fun.value(var));
	    }
	    //TODO funs
	    plotVector(mesh,outputFolder,fileName,v);
	}
	
	public static Vector extendData(Mesh meshFrom, Mesh meshTo, Vector u, double deaultValue) {
		NodeList nodeTo = meshTo.getNodeList();
		int dimTo = nodeTo.size();
		Vector rlt = new SparseVector(dimTo);
		for(int i=1;i<=nodeTo.size();i++) {
			Node node = meshFrom.containNode(nodeTo.at(i));
			if(node != null) {
				rlt.set(i, u.get(node.globalIndex));
			} else {
				rlt.set(i, deaultValue);
			}
		}
		return rlt;
	}	
	/**
	 * 
	 * @param mesh1
	 * @param mesh2
	 * @param u
	 * @return
	 */
	public static Vector extractData(Mesh meshFrom, Mesh meshTo, Vector u) {
		NodeList nodeTo = meshTo.getNodeList();
		int dimTo = nodeTo.size();
		Vector rlt = new SparseVector(dimTo);
		for(int i=1;i<=nodeTo.size();i++) {
			Node node = meshFrom.containNode(nodeTo.at(i));
			if(node != null) {
				rlt.set(i, u.get(node.globalIndex));
			}
		}
		return rlt;
	}	
}
