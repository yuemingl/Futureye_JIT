package p1dssht;

import java.util.HashMap;

import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.algebra.solver.external.SolverJBLAS;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.core.NodeType;
import edu.uta.futureye.core.intf.Assembler;
import edu.uta.futureye.function.MultiVarFunc;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.basic.Vector2Function;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.lib.assembler.AssemblerScalar;
import edu.uta.futureye.lib.element.FELinear1D;
import edu.uta.futureye.lib.weakform.WeakFormLaplace;
import edu.uta.futureye.tutorial.Tools;
import edu.uta.futureye.util.Constant;
import edu.uta.futureye.util.container.ElementList;
import edu.uta.futureye.util.container.NodeList;
import static edu.uta.futureye.function.FMath.*;

public class Driver {
	
	/**  
	 *   This is a 1D steady state heat transfer code for a square geometry using
	 *   the Futureye3.0 FEM library.
	 */
	
	/**
	 * One dimension mesh [0,L]
	 * Point spacing h=L/N
	 * 
	 * @param L: maximum length
	 * @param N: total element number
	 * @return
	 */
	Mesh mesh = null;

	public static Mesh getMesh(double L,double N) {
		Mesh mesh = new Mesh();
		double h=L/N;
		Node node1 = new Node(1,0.0);
		mesh.addNode(node1);
		for(int i=2;i<=N+1;i++) {
			Node node2 = new Node(i,(i-1)*h);
			mesh.addNode(node2);
			NodeList nodeList = new NodeList();
			nodeList.add(node1);
			nodeList.add(node2);
			Element e = new Element(nodeList);
			mesh.addElement(e);
			node1 = node2;
		}
		return mesh;
		
	}
	
	public static void main(String[] args) {
		
		//Generate mesh which just stores nodes and elements
		final double L = .4;
		int N = 40; 
		int maxIter = 60;
		
		Mesh mesh = getMesh(L,N);
		Mesh meshExact = getMesh(L,5*N);
		
        //basic relationship between nodes and elements
        mesh.computeNodeBelongsToElements();

        // Right Dirichlet boundary condition:
        MathFunc BCR = C(400.0);
        
        // Left Neumann boundary condition:
        MathFunc BCL = C0;

        //Right hand side(RHS): q = -b*(4300/(b*x^2 + a + 200) + 21/10) + (8600*b^2*x^2)/(b*x^2 + a + 200)^2;
        //  To force T(x) = 1000-3750 x^2 = a + b x^2
        //   With T(0) = 1000, T(0.4) = 400
        //       Tx(0) = 0
        
        MathFunc XX = x.M(x);
        MathFunc q = C(3750.0).M(C(4300.0).D(C(-3750.0).M(XX).A(C(1020))).A(2.1))
        		.A(C(8600.0).M(14062500.0).M(XX).D(C(-3750.0).M(XX).A(C(1200.0)).M(C(-3750.0).M(XX).A(C(1200.0)))));
        System.out.println(q); //Print your expression
        
        //Mark border type
        HashMap<NodeType, MathFunc> mapNTF = new HashMap<NodeType, MathFunc>();
		//Dirichlet boundary of u
		mapNTF.put(NodeType.Dirichlet, new MultiVarFunc("x") {
//			@Override
//			public double apply(Variable v) {
//				double x = v.get("x");
//				//I think the problem is here, you should use absolute value of (x-0.4)
//				if(Math.abs(x-L) < Constant.meshEps) return 1;
//				else return 0;
//			}

			@Override
			public double apply(double... args) {
				double x = args[0];
				//I think the problem is here, you should use absolute value of (x-0.4)
				if(Math.abs(x-L) < Constant.meshEps) return 1;
				else return 0;
			}
		});
        
       mapNTF.put(NodeType.Neumann, null);
        
        mesh.markBorderNode(mapNTF);
        //This function will write the mesh to a file and there is a number corresponding to each node
   	 	//Inner Node  =5
   		//Dirichlet   =1
   		//Neumann     =2
   		//Robin       =3
   		//Hanging Node=10
   		//Unknown     =20
        mesh.writeNodesInfo("ssht\\mesh.dat");
        

        //Use element library to assign degree of freedom (DOF) to element
        ElementList eList = mesh.getElementList();
        FELinear1D feLT = new FELinear1D();
        for(int i=1;i<=eList.size();i++)
        	feLT.assignTo(eList.at(i));

    	WeakFormLaplace weakForm = new WeakFormLaplace();
    	
    	MathFunc k = new FC(1);
    	
    	for(int i=0; i<maxIter; ++i) { //Iterate maxIter times
    		//The function setParam(Function k,Function c,Function g,Function d)
    		//can be used to set parameters for the model: -div(k*div(u)) + c*u = q
    		weakForm.setF(q);
    		
    		//Pass k into the function.
    		weakForm.setParam(k,new FC(0.0),BCL,new FC(1.0)); //Parameters g and d are for Robin boundary condition d*u +  k*u_n = g, set null for Dirichlet condition
    		
    		Assembler assembler = new AssemblerScalar(mesh, weakForm);
    		assembler.assemble();
    		
    		Matrix stiff = assembler.getStiffnessMatrix();
    		Vector load = assembler.getLoadVector();
    		
    		assembler.imposeDirichletCondition(BCR);
    		
    		SolverJBLAS solver = new SolverJBLAS();
    		Vector T = solver.solveDGESV(stiff, load); //Direct solver
    		
    		//Update parameter k
    		final Vector2Function funU = new Vector2Function(T, mesh, "x");
    		//k depends on T, for example k=2*T
    		k = new MultiVarFunc("x") {
    			public double apply(Variable v) {
    				return 1.05 + 2150.0/(funU.apply(v) + 200.0);
    			}

				@Override
				public double apply(double... args) {
					// TODO Auto-generated method stub
					return 0;
				}
    		};
    		//output intermediate results into the folder "ssht", you can use tecplot to open it
    		Tools.plotFunction(mesh, "ssht", String.format("k_iter%d.dat",i), k);
    		
    		//if (i == maxIter-1){
    			Tools.plotVector(mesh, "ssht", String.format("T_iter%d.dat",i), T);
    		//};
    			
    		//True solution: T(x) = 1000-3740x^2
    		Tools.plotFunction(mesh, "ssht", "T_true.dat", C(1000).S(XX.M(3740)));	
    	};
	}
}
