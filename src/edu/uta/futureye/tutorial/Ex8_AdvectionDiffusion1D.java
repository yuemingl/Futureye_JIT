package edu.uta.futureye.tutorial;

import static edu.uta.futureye.function.FMath.C0;

import java.util.HashMap;

import edu.uta.futureye.algebra.SparseVectorHashMap;
import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.algebra.solver.external.SolverJBLAS;
import edu.uta.futureye.core.DOF;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.core.NodeType;
import edu.uta.futureye.function.SingleVarFunc;
import edu.uta.futureye.function.UserDefFunc;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.lib.assembler.AssembleParam;
import edu.uta.futureye.lib.assembler.BasicAssembler;
import edu.uta.futureye.lib.element.FELinearLine1D;
import edu.uta.futureye.lib.weakform.WeakForm;
import edu.uta.futureye.util.Constant;
import edu.uta.futureye.util.MathEx;
import edu.uta.futureye.util.Utils;
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
public class Ex8_AdvectionDiffusion1D {
	
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
	
	public static Vector solve(Mesh mesh, 
			final double L, final int N, double k, double u) {
		// Mark border types
		HashMap<NodeType, MathFunc> mapNTF = new HashMap<NodeType, MathFunc>();
		mapNTF.put(NodeType.Dirichlet, null);
		mesh.markBorderNode(mapNTF);

		// Weak form definition
		WeakForm wf = new WeakForm(
				new FELinearLine1D(),
				(c, v) -> k * c.diff("x") * v.diff("x") + u * c.diff("x") * v,
				(v)    -> C0
				);
		wf.compile();

		// Assembly and boundary condition(s)
		BasicAssembler assembler = new BasicAssembler(wf);
		assembler.assembleGlobal(mesh);
		Matrix stiff = assembler.getGlobalStiffMatrix();
		Vector load = assembler.getGlobalLoadVector();
		// Boundary condition
		Utils.imposeDirichletCondition(stiff, load, mesh, 
			new SingleVarFunc("diri","x") {
				@Override
				public double apply(double... args) {
					double x = args[0];
					if (Math.abs(x) < Constant.meshEps)
						return 0.0;
					else if (Math.abs(x - L) < Constant.meshEps)
						return 1.0;
					else
						return 0.0;
				}
		});

		// Solve linear system
		SolverJBLAS solver = new SolverJBLAS();
		Vector c = solver.solveDGESV(stiff, load);
		System.out.println("c=");
		for (int i = 1; i <= c.getDim(); i++)
			System.out.println(String.format("%.3f", c.get(i)));

		System.out.println("Grid Peclet number=" + (u * L / N) / (2 * k));

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
		// Mark border types
		HashMap<NodeType, MathFunc> mapNTF = new HashMap<NodeType, MathFunc>();
		mapNTF.put(NodeType.Dirichlet, null);
		mesh.markBorderNode(mapNTF);

		MathFunc upwindCoef = new UserDefFunc() {
			//@Override
			public double apply(AssembleParam ap, double... args) {
				DOF dof = ap.element.getAllNodeDOFList().at(ap.testDOFIdx);
				Node node1 = dof.getNodeOwner();
				int index = ap.element.getLocalIndex(node1);
				Node node2 = null;
				if(index == 1) {
					node2 = ap.element.nodes.at(2);
				} else {
					node2 = ap.element.nodes.at(1);
				}
				double coord1 = node1.coord(1);
				double coord2 = node2.coord(1);
				double upwindWeight = 0.0;
				if((coord2-coord1)*u > 0) {
					upwindWeight = -0.1;
				} else {
					upwindWeight = 0.1;
				}
				return upwindWeight;
			}
		};
		// Weak form definition
		WeakForm wf = new WeakForm(
				new FELinearLine1D(),
				//the convection term is weighted with a modified test function
				(c, v) -> k * c.diff("x") * v.diff("x") + u * c.diff("x") * (v + upwindCoef),
				(v)    -> C0
				);
		wf.compile();

		// Assembly and boundary condition(s)
		BasicAssembler assembler = new BasicAssembler(wf);
		assembler.assembleGlobal(mesh);
		Matrix stiff = assembler.getGlobalStiffMatrix();
		Vector load = assembler.getGlobalLoadVector();
		// Boundary condition
		Utils.imposeDirichletCondition(stiff, load, mesh, 
			new SingleVarFunc("diri","x") {
				@Override
				public double apply(double... args) {
					double x = args[0];
					if (Math.abs(x) < Constant.meshEps)
						return 0.0;
					else if (Math.abs(x - L) < Constant.meshEps)
						return 1.0;
					else
						return 0.0;
				}
		});

		// Solve linear system
		SolverJBLAS solver = new SolverJBLAS();
		Vector c = solver.solveDGESV(stiff, load);
		System.out.println("c=");
		for (int i = 1; i <= c.getDim(); i++)
			System.out.println(String.format("%.3f", c.get(i)));

		System.out.println("Grid Peclet number=" + (u * L / N) / (2 * k));

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
		System.out.println("Exact solution:");
		for(int i=1; i<=cReal.getDim(); i++)
			System.out.println((String.format("%.3f", cReal.get(i))));
		return cReal;
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// 1.Generate mesh which just stores nodes and elements
		double L = 1.0;   // Solution domain: [0, L]
		int N = 10;       // Number of elements: 10 25
		double h = L/N;
		Mesh mesh = getMesh(L, N);
		// Compute the basic relationship between nodes and elements
		mesh.computeNodeBelongsToElements();

		double u = 20.0; //velocity
		Vector ctrue      = exactSolution(mesh, L, N, 1.0, u);

		Vector c1       = solve        (mesh, L, N, 1.0, u);
		Vector c2 = solveUpwind  (mesh, L, N, 1.0, u);

		// Optimal upwind method
		double k = 1.0;
		double alpha = u * h / (2 * k);
		double k_tidle = (u * h / 2) * (MathEx.coth(alpha) - 1 / alpha);
		System.out.println("k_tidle=" + k_tidle);
		Vector c3 = solve(mesh, L, N, k + k_tidle, u);
	}

}
