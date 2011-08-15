package edu.uta.futureye.application;

import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.function.AbstractFunction;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.Vector2Function;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.io.MeshReader;
import edu.uta.futureye.lib.element.FELinearTriangle;
import edu.uta.futureye.util.Constant;
import edu.uta.futureye.util.container.ElementList;
import edu.uta.futureye.util.container.NodeList;

public class NitschesMethodTest {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		String gridFileBig = "prostate_test3_ex.grd";
		String gridFileSmall = "prostate_test3.grd";

		ModelDOT model = new ModelDOT();
		model.setMu_a(2.0, 2.6, 0.3, //(x,y;r)
				0.2, //maxMu_a
				1); //type
		model.setDelta(1.5, 3.5);
	
		MeshReader readerForward = new MeshReader(gridFileBig);
		Mesh meshBig = readerForward.read2DMesh();
		MeshReader readerGCM = new MeshReader(gridFileSmall);
		Mesh meshSmall = readerGCM.read2DMesh();
		
		//Use element library to assign degree of freedom (DOF) to element
		ElementList eList = meshBig.getElementList();
		FELinearTriangle linearTriangle = new FELinearTriangle();
		for(int i=1;i<=eList.size();i++)
			linearTriangle.assignTo(eList.at(i));
		meshBig.computeNodeBelongsToElements();
		meshBig.computeNeighborNodes();
		
		eList = meshSmall.getElementList();
		for(int i=1;i<=eList.size();i++)
			linearTriangle.assignTo(eList.at(i));
		meshSmall.computeNodeBelongsToElements();
		meshSmall.computeNeighborNodes();		
		
		Vector uBig = model.solveNeumann(meshBig);
		Tools.plotVector(meshBig, "NitschesMethodTest", "u_big.dat", uBig);

		Vector uSmallForBoundary = Tools.extractData(meshBig, meshSmall, uBig);
		NodeList nodes = meshSmall.getNodeList();
		//NOT necessary:
		for(int i=1;i<=nodes.size();i++) {
			if(nodes.at(i).isInnerNode())
				uSmallForBoundary.set(i,0.0);
		}
		Vector uSmall = model.solveDirichlet(meshSmall, new Vector2Function(uSmallForBoundary));
		Tools.plotVector(meshSmall, "NitschesMethodTest", "u_small.dat", uSmall);
		
		
		Vector uNitsches = model.solveNitsches(meshSmall, new Vector2Function(uSmallForBoundary), 0.00001);
		Tools.plotVector(meshSmall, "NitschesMethodTest", "u_Nitsches.dat", uNitsches);
		
		
		
	}

}
