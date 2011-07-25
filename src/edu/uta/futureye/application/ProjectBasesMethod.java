package edu.uta.futureye.application;

import java.io.File;
import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;

import edu.uta.futureye.algebra.CompressedColMatrix;
import edu.uta.futureye.algebra.CompressedRowMatrix;
import edu.uta.futureye.algebra.FullVector;
import edu.uta.futureye.algebra.Solver;
import edu.uta.futureye.algebra.SparseMatrix;
import edu.uta.futureye.algebra.SparseVector;
import edu.uta.futureye.algebra.intf.AlgebraMatrix;
import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.core.NodeType;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.io.MeshReader;
import edu.uta.futureye.io.MeshWriter;
import edu.uta.futureye.lib.element.FELinearTriangle;
import edu.uta.futureye.util.Constant;
import edu.uta.futureye.util.PairDoubleInteger;
import edu.uta.futureye.util.container.ElementList;
import edu.uta.futureye.util.container.NodeList;

public class ProjectBasesMethod {
	protected static String outputFolder = "ProjectBasesMethod";
	public boolean debug = false;
	Mesh meshBig;
	Mesh mesh;
	
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
	
	public void run() {
		readMeshTriangle();
		NodeList upSideNodes = getUpsideNodes(mesh);
		
		ModelDOT model = new ModelDOT();
		model.setDelta(2.0, 3.5);
		
		//Generate bases on upside
		double h = 0.2;
		List<Vector> uSmalls = new ArrayList<Vector>();
		for(int xi=1;xi<=100;xi++) {
			double x = 1.0+(xi-1)*h;
			if(x>4.0) break;
			model.setMu_a(x, 2.8, 0.3, 1.0, 1);
			plotFunction(meshBig, model.mu_a, String.format("mu_a%02d.dat",xi));
			Vector u = model.solveNeumann(meshBig);
			plotVector(meshBig, u, String.format("u_big%02d.dat",xi));
			//截取meshBig的部分解到mesh上
			Vector uSmall = Tools.extractData(meshBig, mesh, u);
			plotVector(mesh, uSmall, String.format("u_%02d.dat",xi));
			uSmalls.add(uSmall);
		}
		
		//Generate measurement data on upside
		model.setMu_a(3.1, 2.8, 0.3, 1.0, 1);
		plotFunction(meshBig, model.mu_a, String.format("mu_a.dat"));
		Vector ua = model.solveNeumann(meshBig);
		plotVector(meshBig, ua, String.format("u_big.dat"));
		//截取meshBig的部分解到mesh上
		Vector uaSamll = Tools.extractData(meshBig, mesh, ua);
		plotVector(mesh, uaSamll, String.format("u.dat"));
		
		int nRow = upSideNodes.size();
		int nCol = uSmalls.size();
		
		SparseMatrix A = new SparseMatrix(nRow,nCol);
		SparseVector f = new SparseVector(nRow);
		for(int c=1;c<=nCol;c++) {
			Vector base = uSmalls.get(c-1);
			for(int r=1;r<=nRow;r++) {
				A.set(r, c, 
						base.get(upSideNodes.at(r).globalIndex));
			}
		}
		for(int r=1;r<=nRow;r++) {
			f.set(r, uaSamll.get(upSideNodes.at(r).globalIndex));
		}
		
		//A.print();
		//f.print();
		
		SparseVector x0 = new SparseVector(nCol);
		x0.set(10, 1.0);
		x0.set(11, 1.0);
		x0.set(12, 1.0);
		
		CompressedColMatrix AA = new CompressedColMatrix(A, false);
		CompressedRowMatrix C = new CompressedRowMatrix(); //
		AlgebraMatrix AAT = AA.getTrans();
		AAT.mult(AA, C);
		FullVector ff = new FullVector(f);
		FullVector g  = new FullVector(nCol);
		AAT.mult(ff, g);
		
		SparseMatrix diagLmd = new SparseMatrix(nRow,nRow);
		for(int i=1;i<=nRow;i++)
			diagLmd.set(i, i, 1.5/4.0);
		C.axpy(1.0, new CompressedRowMatrix(diagLmd, false));
		g.add(new FullVector(x0));
		
		FullVector x  = new FullVector(x0);
		Solver sol = new Solver();
		sol.epsIter = 1e-10;
		sol.solveCGS(C, g, x);
		//C.print();
		//g.print();
		x.print();

	}
	
	
	public void run2() {
		readMeshTriangle();
		NodeList upSideNodes = getUpsideNodes(mesh);
		
		ModelDOT model = new ModelDOT();
		model.setDelta(2.0, 3.5);
		
		//Generate bases on upside
		List<Vector> uSmalls = new ArrayList<Vector>();
		double []xx = {2.8,2.9,3.0,3.1,3.2};
		double []yy = {2.6,2.7,2.8,2.9,3.0};
		for(int i=0;i<xx.length;i++) {
			for(int j=0;j<yy.length;j++) {
				int cnt = i*yy.length + j;
				model.setMu_a(xx[i], yy[j], 0.2, 1.0, 1);
				plotFunction(meshBig, model.mu_a, String.format("mu_a%02d.dat",cnt));
				Vector u = model.solveNeumann(meshBig);
				plotVector(meshBig, u, String.format("u_big%02d.dat",cnt));
				//截取meshBig的部分解到mesh上
				Vector uSmall = Tools.extractData(meshBig, mesh, u);
				plotVector(mesh, uSmall, String.format("u_%02d.dat",cnt));
				uSmalls.add(uSmall);
			}
		}
		
		//Generate measurement data on upside
		model.setMu_a(2.83, 2.63, 0.3, 1.0, 1);
		plotFunction(meshBig, model.mu_a, String.format("mu_a.dat"));
		Vector ua = model.solveNeumann(meshBig);
		plotVector(meshBig, ua, String.format("u_big.dat"));
		//截取meshBig的部分解到mesh上
		Vector uaSamll = Tools.extractData(meshBig, mesh, ua);
		plotVector(mesh, uaSamll, String.format("u.dat"));
		
		int nRow = upSideNodes.size();
		int nCol = uSmalls.size();
		
		SparseMatrix A = new SparseMatrix(nRow,nCol);
		SparseVector f = new SparseVector(nRow);
		for(int c=1;c<=nCol;c++) {
			Vector base = uSmalls.get(c-1);
			for(int r=1;r<=nRow;r++) {
				A.set(r, c, 
						base.get(upSideNodes.at(r).globalIndex));
			}
		}
		for(int r=1;r<=nRow;r++) {
			f.set(r, uaSamll.get(upSideNodes.at(r).globalIndex));
		}
		
		//A.print();
		//f.print();
		
		SparseVector x0 = new SparseVector(nCol);
		x0.set(1, 1.0);
		
		CompressedColMatrix AA = new CompressedColMatrix(A, false);
		CompressedRowMatrix C = new CompressedRowMatrix(); //
		AlgebraMatrix AAT = AA.getTrans();
		AAT.mult(AA, C);
		FullVector ff = new FullVector(f);
		FullVector g  = new FullVector(nCol);
		AAT.mult(ff, g);
		
		SparseMatrix diagLmd = new SparseMatrix(nRow,nRow);
		for(int i=1;i<=nRow;i++)
			diagLmd.set(i, i,10.0/4.0);
		C.axpy(1.0, new CompressedRowMatrix(diagLmd, false));
		g.add(new FullVector(x0));
		
		FullVector x  = new FullVector(x0);
		Solver sol = new Solver();
		sol.epsIter = 1e-10;
		sol.solveCGS(C, g, x);
		//C.print();
		//g.print();
		x.print();
		
		
		Function rlt_mu_a = FC.c0;
		double[] coef = x.getData();
		for(int i=0;i<xx.length;i++) {
			for(int j=0;j<yy.length;j++) {
				int cnt = i*yy.length + j;
				model.setMu_a(xx[i], yy[j], 0.2, 1.0, 1);
				rlt_mu_a = rlt_mu_a.A(model.mu_a.M(coef[cnt]));
			}
		}
		plotFunction(meshBig, rlt_mu_a, String.format("mu_a_rlt.dat"));

	}
	
	public NodeList getUpsideNodes(Mesh mesh) {
		NodeList upSideNodes = new NodeList();
		NodeList nodes = mesh.getNodeList();
		for(int i=1;i<=nodes.size();i++) {
			Node node = nodes.at(i);
			double y = node.coord(2);
			if(Math.abs(y-3.0)<Constant.meshEps) {
				upSideNodes.add(node);
			}
		}
		List<Node> list = upSideNodes.toList();
		Collections.sort(list, new Comparator<Node>() {
			@Override
			public int compare(Node o1, Node o2) {
				return o1.coord(1) > o2.coord(1) ? 1 : -1;
		}});
		return upSideNodes;
	}
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		ProjectBasesMethod pbm = new ProjectBasesMethod();
		pbm.run2();
	}

}
