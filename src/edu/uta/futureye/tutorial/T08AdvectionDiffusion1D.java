package edu.uta.futureye.tutorial;

import java.util.HashMap;

import edu.uta.futureye.algebra.SparseVectorHashMap;
import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.algebra.solver.external.SolverJBLAS;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.core.NodeType;
import edu.uta.futureye.core.intf.WeakForm;
import edu.uta.futureye.function.AbstractMathFun;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.function.intf.ScalarShapeFunction;
import edu.uta.futureye.io.MeshWriter;
import edu.uta.futureye.lib.assembler.AssemblerScalar;
import edu.uta.futureye.lib.element.FELinear1D;
import edu.uta.futureye.lib.weakform.WeakFormAdvectionDiffusion1D;
import edu.uta.futureye.lib.weakform.WeakFormBuilder;
import edu.uta.futureye.util.Constant;
import edu.uta.futureye.util.MathEx;
import edu.uta.futureye.util.container.ElementList;
import edu.uta.futureye.util.container.NodeList;

/**
 * <blockquote><pre>
 * One-dimensional advection-diffusion problem:
 * 
 *   u*c_x = k*c_xx
 *   c=0 at x=0
 *   c=1 at x=L
 *   
 * where
 *   x in [0,L]
 *   c=c(x): particles or energy(e.g. salt density, Heat...) are transferred inside 
 *                 a physical system due to two processes: diffusion and convection
 *   u: flow velocity
 *   k: diffusivity
 * 
 * Weak form:
 *   (k*c_x,v_x) + (u*c_x,v) = 0, for all v
 * 
 * Real solution:
 * 
 *   c=(1-e^(Pe*x/L))/(1-e^Pe)
 * where
 *   Pe=u*L/k (gobal Peclet number)
 * 
 * </blockquote></pre>
 * 
 * @author liuyueming
 *
 */
public class T08AdvectionDiffusion1D {
	
	/**
	 * One dimension mesh [0,L]
	 * Point spacing h=L/N
	 * 
	 * @param L: maximum length
	 * @param N: total element number
	 * @return
	 */
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
	
	/**
	 * 
	 * @param mesh
	 * @param L
	 * @param N
	 * @param k diffusivity
	 * @param u flow velocity
	 * @return
	 */
	public static Vector solve(Mesh mesh, 
			final double L, final int N, double k, double u) {
        //2.Mark border types
        HashMap<NodeType, MathFunc> mapNTF =
                new HashMap<NodeType, MathFunc>();
        mapNTF.put(NodeType.Dirichlet, null);
        mesh.markBorderNode(mapNTF);

        //3.Use element library to assign degrees of
        //  freedom (DOF) to element
        ElementList eList = mesh.getElementList();
        FELinear1D feLT = new FELinear1D();
        for(int i=1;i<=eList.size();i++)
            feLT.assignTo(eList.at(i));

        //4.Weak form
		WeakFormBuilder wfb = new WeakFormBuilder() {
			@Override
			public MathFunc makeExpression(Element e, Type type) {
				ScalarShapeFunction u = getScalarTrial();
				ScalarShapeFunction v = getScalarTest();
				MathFunc fk = getParam("k",e);
				MathFunc fu = getParam("u",e);
				switch(type) {
					case LHS_Domain:
						return fk.M(u.diff("x").M(v.diff("x"))).A(fu.M(u.diff("x").M(v)));
					default:
						return null;
				}
			}
		};
		wfb.addParam("k", FC.c(k));
		wfb.addParam("u", FC.c(u));
		WeakForm wf = wfb.getScalarWeakForm();
        
        //5.Assembly process
        AssemblerScalar assembler =
                new AssemblerScalar(mesh, wf);
        assembler.assemble();
        Matrix stiff = assembler.getStiffnessMatrix();
        Vector load = assembler.getLoadVector();
        //Boundary condition
        assembler.imposeDirichletCondition(new AbstractMathFun("x") {
        	@Override
        	public double apply(Variable v) {
        		double x = v.get("x");
        		if(Math.abs(x)<Constant.meshEps)
        			return 0.0;
        		else if(Math.abs(x-L)<Constant.meshEps)
        			return 1.0;
        		else
        			return 0.0;
        	}
        });

        //6.Solve linear system
        SolverJBLAS solver = new SolverJBLAS();
        Vector c = solver.solveDGESV(stiff, load);
//        System.out.println("c=");
//        for(int i=1;i<=c.getDim();i++)
//            System.out.println(String.format("%.3f", c.get(i)));
        System.out.println("Grid Peclet number="+(u*L/N)/(2*k));
        
        return c;
	}
	
	/**
	 * 
	 * @param mesh
	 * @param L
	 * @param N
	 * @param k diffusivity
	 * @param u flow velocity
	 * @return
	 */
	public static Vector solveUpwind(Mesh mesh, 
			final double L, final int N, double k, double u) {
        //2.Mark border types
        HashMap<NodeType, MathFunc> mapNTF =
                new HashMap<NodeType, MathFunc>();
        mapNTF.put(NodeType.Dirichlet, null);
        mesh.markBorderNode(mapNTF);

        //3.Use element library to assign degrees of
        //  freedom (DOF) to element
        ElementList eList = mesh.getElementList();
        FELinear1D feLT = new FELinear1D();
        for(int i=1;i<=eList.size();i++)
            feLT.assignTo(eList.at(i));

        //4.Weak form
        WeakFormAdvectionDiffusion1D wf = new WeakFormAdvectionDiffusion1D();
        wf.setParam(FC.c(k));
        wf.setConvectionVelocity(FC.c(u));
        wf.setF(FC.C0);
        //5.Assembly process
        AssemblerScalar assembler =
                new AssemblerScalar(mesh, wf);
        assembler.assemble();
        Matrix stiff = assembler.getStiffnessMatrix();
        Vector load = assembler.getLoadVector();
        //Boundary condition
        assembler.imposeDirichletCondition(new AbstractMathFun("x") {
        	@Override
        	public double apply(Variable v) {
        		double x = v.get("x");
        		if(Math.abs(x)<Constant.meshEps)
        			return 0.0;
        		else if(Math.abs(x-L)<Constant.meshEps)
        			return 1.0;
        		else
        			return 0.0;
        	}
        });

        //6.Solve linear system
        SolverJBLAS solver = new SolverJBLAS();
        Vector c = solver.solveDGESV(stiff, load);
//        System.out.println("c=");
//        for(int i=1;i<=c.getDim();i++)
//            System.out.println(String.format("%.3f", c.get(i)));
        System.out.println("Grid Peclet number="+(u*L/N)/(2*k));
        
        return c;
	}
	
	public static Vector exactSolution(Mesh mesh, 
			final double L, final int N, double k, double u) {
		Vector cReal = new SparseVectorHashMap(N+1);
		double Pe = u*L/k;
		for(int i=1;i<=N+1;i++) {
			double x = mesh.getNodeList().at(i).coord(1);
			cReal.set(i, (1.0-Math.exp(Pe*x/L))/(1-Math.exp(Pe)));
		}
		return cReal;
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		//1.Generate mesh which just stores nodes and elements
		double L = 1.0;
		int N = 10; //=10 25
		Mesh mesh = getMesh(L,N);
		Mesh meshExact = getMesh(L,5*N);
        //basic relationship between nodes and elements
        mesh.computeNodeBelongsToElements();
        
		
        Vector c1 = solve(mesh,L,N,1.0,10);
		Vector c2 = solve(mesh,L,N,1.0,20);
		Vector c2upwind = solveUpwind(mesh,L,N,1.0,20);
		Vector c3 = solve(mesh,L,N,1.0,50);
		Vector c3upwind = solveUpwind(mesh,L,N,1.0,50);
		
		
		Vector ce1 = exactSolution(meshExact,L,5*N,1.0,10);
		Vector ce2 = exactSolution(meshExact,L,5*N,1.0,20);
		Vector ce3 = exactSolution(meshExact,L,5*N,1.0,50);
		
		//Optimal upwind method
        double k=1.0;
        double u=50;
		double h = L/N;
		double alpha = u*h/(2*k);
		double k_tidle = (u*h/2)*(MathEx.coth(alpha)-1/alpha);
		System.out.println("k_tidle="+k_tidle);
        Vector c31 = solve(mesh,L,N,k+k_tidle,u);
        
        //7.Output results to an Techplot format file
        MeshWriter writer = new MeshWriter(mesh);
        MeshWriter writerEx = new MeshWriter(meshExact);
        writer.writeTechplot("./tutorial/AdvectionDiffusion1D.dat", c1,c2,c2upwind,c3,c3upwind,c31);
        writerEx.writeTechplot("./tutorial/AdvectionDiffusion1Dexact.dat", ce1,ce2,ce3);
	}

}
