package edu.uta.futureye.test;

import java.util.HashMap;

import edu.uta.futureye.algebra.SchurComplementSolver;
import edu.uta.futureye.algebra.SparseVector;
import edu.uta.futureye.algebra.intf.BlockMatrix;
import edu.uta.futureye.algebra.intf.BlockVector;
import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.core.DOF;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.NodeType;
import edu.uta.futureye.function.basic.FAxpb;
import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.basic.FX;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.io.MeshReader;
import edu.uta.futureye.io.MeshWriter;
import edu.uta.futureye.lib.assembler.AssemblerMixedLaplace;
import edu.uta.futureye.lib.shapefun.RaviartThomas2D0;
import edu.uta.futureye.lib.weakform.WeakFormMixedLaplace;
import edu.uta.futureye.util.container.EdgeList;
import edu.uta.futureye.util.container.ElementList;
import edu.uta.futureye.util.container.NodeList;

public class MixedLaplaceTest {
	public static void main(String[] args) {
		//String meshName = "patch_EBmfem";//OK
		//String meshName = "patch_triangle";//OK
		//String meshName = "patch_triangle2";//OK
		//String meshName = "patch_triangle3";//OK
		//String meshName = "patch_triangle4";//NO
		//String meshName = "triangle"; //NO
		String meshName = "triangle2"; //OK
		
		MeshReader reader = new MeshReader(meshName+".grd");
		Mesh mesh = reader.read2DMesh();
		
		mesh.computeNodeBelongsToElements();
		mesh.computeGlobalEdge();
				EdgeList edgeList = mesh.getEdgeList();
		System.out.println(edgeList.size());
		for(int i=1;i<=edgeList.size();i++) {
			System.out.println(edgeList.at(i));
		}
	    ElementList eList = mesh.getElementList();
//	    for(int i=1;i<=eList.size();i++) {
//	    	EdgeList<EdgeLocal> egList = eList.at(i).edges;
//	    	for(int j=1;j<=egList.size();j++) {
//	    		EdgeLocal eg = egList.at(j);
//	    		double sign = eg.getNormVector().dot(eg.getGlobalEdge().getNormVector());
//	    		System.out.println(sign + "  " + eg.toString() + " <--> " + eg.getGlobalEdge());
//	    	}
//	    }
		
		
		HashMap<NodeType, Function> mapNTF = new HashMap<NodeType, Function>();
		mapNTF.put(NodeType.Dirichlet, null);
		mesh.markBorderNode(mapNTF);

		RaviartThomas2D0[] shapeFun = new RaviartThomas2D0[3];
		for(int i=0;i<3;i++)
			shapeFun[i] = new RaviartThomas2D0(i+1);
		
		//Asign degree of freedom to element
		ElementList eleList = mesh.getElementList();
		int nElements = eleList.size();
		int nEdges = mesh.getEdgeList().size();
		//int[] order = {0,1,3,2};
		for(int i=1;i<=nElements;i++) {
			Element e = eleList.at(i);
			int nEleEdges = e.edges().size();
			for(int j=1;j<=nEleEdges;j++) {
				//Asign shape function to DOF
				DOF dof = new DOF(
						j,//order[j],
						e.edges().at(j).getGlobalEdge().getGlobalIndex(),
						shapeFun[j-1]);
				//Associate DOF to edge
				e.addEdgeDOF(j, dof);
			}
			e.addVolumeDOF(new DOF(
					1,
					nEdges + i,
					null));
		}
		
		WeakFormMixedLaplace weakForm = new WeakFormMixedLaplace();
		
		//-\Delta{u} = f
		//u(x,y)=0, (x,y)\in\partial{\Omega}
		Function fxm5 = new FAxpb("x",1.0,-5.0);
		Function fym5 = new FAxpb("y",1.0,-5.0);
		weakForm.setF(
				//new FConstant(1.0)
	
				//\Omega = [0,10]*[0,10]
				//u=[(x-5)^2-25]*[(y-5)^2-25]
				//f=-2*( (x-5)^2 + (y-5)^2 ) + 100
//				FOBasic.Plus(
//					FOBasic.Plus(
//						FOBasic.Mult(new FConstant(-2.0), FOBasic.Power(fxm5, new FConstant(2.0)) ),
//						FOBasic.Mult(new FConstant(-2.0), FOBasic.Power(fym5, new FConstant(2.0)) )
//						),new FConstant(100.0)
//					)
				
				//\Omega = [-3,3]*[-3,3]
				//u=(x^2-9)*(y^2-9)
				//f=-2*(x^2+y^2)+36
				FC.c(-2.0).M(
						FX.fx.M(FX.fx).A(FX.fy.M(FX.fy))
					).A(FC.c(36.0))
					
				);
		weakForm.setParam(
					null,
					null,
					FC.c(6.0).M(FX.fy.M(FX.fy)).S(FC.c(54.0)),
					null //Robin: 6*y^2-54
				);
		
		AssemblerMixedLaplace assembler = new AssemblerMixedLaplace(mesh, weakForm);
		System.out.println("Begin Assemble...");
		long begin = System.currentTimeMillis();
		assembler.assemble();
		Matrix stiff = assembler.getStiffnessMatrix();
		Vector load = assembler.getLoadVector();
		//assembler.imposeDirichletCondition(new FC(0.0));
		long end = System.currentTimeMillis();
		System.out.println("Assemble done!");
		System.out.println("Time used:"+(end-begin));
		
		SchurComplementSolver solver = new SchurComplementSolver(
				(BlockMatrix)stiff, 
				(BlockVector)load);
		//stiff.print();
		//load.print();
		Vector u = solver.solve();
	    System.out.println("u=");
	    for(int i=1;i<=u.getDim();i++)
	        System.out.println(String.format("%.4f", u.get(i)));
	    
	    //u=(Flux Disp)
	    Vector disp = ((BlockVector)u).getBlock(2);
	    eList = mesh.getElementList();
	    Vector out_disp = new SparseVector(mesh.getNodeList().size());
	    for(int i=1;i<=eList.size();i++) {
	    	NodeList nList = eList.at(i).nodes;
	    	for(int j=1;j<=nList.size();j++) {
	    		if(disp.get(eList.at(i).globalIndex) > out_disp.get(nList.at(j).globalIndex))
	    		out_disp.set(
	    				nList.at(j).globalIndex,
	    				disp.get(eList.at(i).globalIndex)
	    				);
	    	}
	    }
	    
	    MeshWriter writer = new MeshWriter(mesh);
	    writer.writeTechplot(meshName+"_mixed_out.dat", out_disp);

	}
}
