package edu.uta.futureye.tutorial;

import java.util.HashMap;

import edu.uta.futureye.algebra.Solver;
import edu.uta.futureye.algebra.SparseMatrix;
import edu.uta.futureye.algebra.intf.BlockVector;
import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.NodeType;
import edu.uta.futureye.function.AbstractFunction;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.basic.SpaceVectorFunction;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.io.MeshReader;
import edu.uta.futureye.io.MeshWriter;
import edu.uta.futureye.lib.assembler.AssemblerVector;
import edu.uta.futureye.lib.element.FELinearTriangleVector;
import edu.uta.futureye.lib.weakform.WeakFormElasticIsoPlaneStress2D;
import edu.uta.futureye.util.Constant;
import edu.uta.futureye.util.container.ElementList;
import edu.uta.futureye.util.container.NodeList;
import edu.uta.futureye.util.container.ObjIndex;

public class PlaneElasticDam {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
      MeshReader reader = new MeshReader("elastic_dam.grd");
      Mesh mesh = reader.read2DMesh();
      //Compute geometry relationship of nodes and elements
      mesh.computeNodeBelongsToElements();

      //2.Mark border types
      HashMap<NodeType, Function> mapNTF =
              new HashMap<NodeType, Function>();
      mapNTF.put(NodeType.Robin, new AbstractFunction("x","y"){
      	@Override
      	public double value(Variable v) {
      		double x = v.get("x");
      		double y = v.get("y");
      		//y=-2x+20
      		if(Math.abs(y-(-2.0*x+20))<Constant.meshEps)
      			return 1;
      		else
      			return 0;
      	}
      });
      mapNTF.put(NodeType.Dirichlet, new AbstractFunction("x","y"){
      	@Override
      	public double value(Variable v) {
      		double x = v.get("x");
      		if(Math.abs(x)<Constant.meshEps)
      			return 1;
      		else
      			return 0;
      	}
      });
      //注意：位移向量的每个分量都需要制定边界条件
      mesh.markBorderNode(new ObjIndex(1,2),mapNTF);

      NodeList nodes = mesh.getNodeList();
      ElementList eles = mesh.getElementList();
      for(int i=1;i<=nodes.size();i++) {
      	if(nodes.at(i).getNodeType()==NodeType.Robin)
      		System.out.println(NodeType.Robin+":"+nodes.at(i));
      	if(nodes.at(i).getNodeType()==NodeType.Dirichlet)
      		System.out.println(NodeType.Dirichlet+":"+nodes.at(i));
     	
      }
      for(int i=1;i<=eles.size();i++) {
      	System.out.println(eles.at(i));
      }

      //3.Use element library to assign degrees of
      //  freedom (DOF) to element
      ElementList eList = mesh.getElementList();
//      FEBilinearRectangleVector fe = new FEBilinearRectangleVector();
      FELinearTriangleVector fe = new FELinearTriangleVector();
      fe.initDOFIndexGenerator(mesh.getNodeList().size());
      for(int i=1;i<=eList.size();i++)
          fe.assignTo(eList.at(i));

      //4.Weak form
      WeakFormElasticIsoPlaneStress2D weakForm = new WeakFormElasticIsoPlaneStress2D();
      SpaceVectorFunction b = new SpaceVectorFunction(2);
      SpaceVectorFunction t = new SpaceVectorFunction(2);
      b.set(1,FC.c0);
      b.set(2,FC.c0);
      t.set(1,FC.c0);
      t.set(2,new AbstractFunction("x","y"){
        	@Override
          	public double value(Variable v) {
          		double y = v.get("y");
          		return -(10-y);
          	}
          });
      weakForm.setF(b,t);

      //5.Assembly process
      AssemblerVector assembler =
              new AssemblerVector(mesh, weakForm, fe);
      System.out.println("Begin Assemble...");
      assembler.assemble();
      Matrix stiff = assembler.getStiffnessMatrix();
      Vector load = assembler.getLoadVector();
      //Boundary condition
      SpaceVectorFunction diri = new SpaceVectorFunction(2);
      diri.set(1,FC.c0);
      diri.set(2,FC.c0);
      assembler.imposeDirichletCondition(diri);
      System.out.println("Assemble done!");

      //6.Solve linear system
      Solver solver = new Solver();
      SparseMatrix stiff2 = new SparseMatrix(stiff.getRowDim(),stiff.getColDim());
      stiff2.setAll(0, 0, stiff.getAll());
      Vector u = solver.solveCGS(stiff2, load);
      System.out.println("u=");
      for(int i=1;i<=u.getDim();i++)
          System.out.println(String.format("%.3f", u.get(i)));

      BlockVector blkU = (BlockVector)u;
      
      //7.Output results to an Techplot format file
      MeshWriter writer = new MeshWriter(mesh);
      writer.writeTechplot("tuitoral_ElasticDam.dat", blkU.getBlock(1),
      		blkU.getBlock(2));	}

}
