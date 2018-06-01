package edu.uta.futureye.tutorial;

import static edu.uta.futureye.function.FMath.C;
import static edu.uta.futureye.function.FMath.C0;
import static edu.uta.futureye.function.FMath.grad;

import java.util.HashMap;

import edu.uta.futureye.algebra.SparseBlockMatrix;
import edu.uta.futureye.algebra.SparseBlockVector;
import edu.uta.futureye.algebra.SparseMatrixRowMajor;
import edu.uta.futureye.algebra.SparseVectorHashMap;
import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.SparseVector;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.algebra.solver.SchurComplementStokesSolver;
import edu.uta.futureye.core.Edge;
import edu.uta.futureye.core.EdgeLocal;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.core.NodeLocal;
import edu.uta.futureye.core.NodeType;
import edu.uta.futureye.core.Vertex;
import edu.uta.futureye.core.intf.FiniteElement;
import edu.uta.futureye.core.intf.VecFiniteElement;
import edu.uta.futureye.function.MultiVarFunc;
import edu.uta.futureye.function.UserDefFunc;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.SpaceVectorFunction;
import edu.uta.futureye.function.basic.Vector2MathFunc;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.function.intf.VecMathFunc;
import edu.uta.futureye.io.MeshReader;
import edu.uta.futureye.lib.assembler.AssembleParam;
import edu.uta.futureye.lib.assembler.BasicVecAssembler;
import edu.uta.futureye.lib.element.FEQuadraticV_LinearP;
import edu.uta.futureye.lib.weakform.VecWeakForm;
import edu.uta.futureye.util.Constant;
import edu.uta.futureye.util.FutureyeException;
import edu.uta.futureye.util.Utils;
import edu.uta.futureye.util.container.ElementList;
import edu.uta.futureye.util.container.ObjIndex;
import edu.uta.futureye.util.container.ObjList;


/**
 * Problem: 2D Navier-Stokes, lid-driven cavity
 *
 * Problem:
 * -\nabla({k\nabla{\vec{u}}})
 * 		+ \vec{U}\cdot\nabla\vec{u}
 * 		+ c\vec{u} 
 * 		+ \nabla{p} 
 * 		= \vec{f}
 * div{\vec{u}} = 0
 * 
 * Weak form:
 *   find \vec{u} \in H_0^1(div;\Omega), p \in L_2(\Omega)
 *   such that, for all \vec{v} \in H_0^1(div;\Omega), q \in L_2(\Omega)
 *   
 *   (\nabla\vec{v},k*\nabla\vec{u})
 *   		+ (\vec{U}\cdot\nabla\vec{u},\vec{v})
 *   		+ (c*\vec{u},\vec{v})
 *   		- (div{\vec{v}},p)
 *  		+ (q,div{\vec{u}}) 
 *			= (\vec{v},\vec{f})
 *
 *   k* [(v1_x,u1_x) + (v1_y,u1_y) + (v2_x,u2_x) + (v2_y,u2_y) ]
 *   		+ [(U1*u1_x,v1)+(U2*u1_y,v1)] + [(U1*u2_x,v2)+(U2*u2_y,v2)]
 *   		+ c*[(u1,v1)+(u2,v2)]
 *			- (v1_x+v2_y,p)
 *			+ (q,u1_x+u2_y)
 *			= (v1,f1)+(v2,f2)
 *
 * where
 *   \vec{u}=(u1,u2): velocity vector field
 *   \vec{f}=(f1,f2): body force
 *   \vec{U}=(U1,U2): previous velocity
 * </blockquote></pre>
 * 
 * Ref.
 * 1. T.E. Teaduyar Stabilized Finite Element Formulations for Incompressible Flow Computations
 * 
 * 2. Rajiv Sampath and Nicholas Zabaras Design and object-oriented implementation of a 
 *    preconditioned-stabilized incompressible NaiverStokes solver using equal-order-interpolation 
 *    velocity pressure elements
 *
 */
public class Ex11_NavierStokesBox {
	protected String outputFolder = "NavierStokesBox";
	protected String file = null;
	
	protected Mesh mesh = null;
	protected Mesh meshOld = null;
	
	//Dirichlet boundary condition
	protected VecMathFunc diri = null;
	
	//viscosity
	protected double nu = 0.001; 
	
	int maxNonlinearIter = 30;
	double nonlinearError = 1e-2;
	
	/**
	 * 
	 * @param testCaseNo
	 */
	public void init(int testCaseNo) {
		//Read mesh from input file
		if(testCaseNo == 1)
			file = "stokes_box";//[-3,3]*[-3,3]
		else if(testCaseNo == 2) 
			file = "stokes_box2";//[0,1]*[0,1] 64*64
		else if(testCaseNo == 3) 
			file = "stokes_box3";//[0,1]*[0,1] 32*32
		else
			throw new FutureyeException("testCaseNo should be 1,2");
		
		MeshReader reader = new MeshReader("grids/"+file+".grd");
		MeshReader reader2 = new MeshReader("grids/"+file+".grd");
		mesh = reader.read2DMesh();
		meshOld = reader2.read2DMesh();
		mesh.nVertex = mesh.getNodeList().size();
		
		//Add nodes for quadratic element, original mesh stored in meshOld
		if(testCaseNo == 1) {
			for(int i=1;i<=mesh.getElementList().size();i++) {
				Element e = mesh.getElementList().at(i);
				e.adjustVerticeToCounterClockwise();
				ObjList<EdgeLocal> edges = e.edges();
				int nNode = e.nodes.size();
				for(int j=1;j<=edges.size();j++) {
					EdgeLocal edge = edges.at(j);
					Vertex l = edge.beginVertex();
					Vertex r = edge.endVertex();
					double cx = (l.coord(1)+r.coord(1))/2.0;
					double cy = (l.coord(2)+r.coord(2))/2.0;
					Node node = new Node(mesh.getNodeList().size()+1, cx,cy);
					Node findNode = mesh.findNode(node);
					if(findNode == null) {
						edge.addEdgeNode(new NodeLocal(++nNode,node));
						mesh.addNode(node);
					} else {
						edge.addEdgeNode(new NodeLocal(++nNode,findNode));
					}
				}
				e.applyChange();
			}
		}
		//Geometry relationship
		mesh.computeNodeBelongsToElements();
		
		//Mark border type
		HashMap<NodeType, MathFunc> mapNTF_uv = new HashMap<NodeType, MathFunc>();
		mapNTF_uv.put(NodeType.Dirichlet, null);
		
		HashMap<NodeType, MathFunc> mapNTF_p = new HashMap<NodeType, MathFunc>();
		mapNTF_p.put(NodeType.Neumann, null);
		
		mesh.markBorderNode(new ObjIndex(1,2),mapNTF_uv);
		mesh.markBorderNode(3,mapNTF_p);

		diri = new SpaceVectorFunction(3);
		if(testCaseNo == 1)
			diri.set(1, new MultiVarFunc("diri", "x", "y") {
						@Override
						public double apply(double... args) {
							double y = args[this.argIdx[1]];
							if(Math.abs(y-3.0)<Constant.meshEps)
								return 1.0;
							else
								return 0.0;
						}
					});
		else
			diri.set(1, new MultiVarFunc("diri", "x", "y") {
				@Override
				public double apply(double... args) {
					double y = args[this.argIdx[1]];
					if(Math.abs(y-1.0)<Constant.meshEps)
						return 1.0;
					else
						return 0.0;
				}
			});
		diri.set(2, C0);
		diri.set(3, C0);
	}

	/**
	 * Solve the steady Navier-Stokes equation iteratively 
	 * by passing the previous solution uk as the convection speed.
	 * 
	 * @param nIter
	 * @param uk
	 * @return
	 */
	public SparseBlockVector nonlinearIterSteady(int nIter, SpaceVectorFunction uk) {
		// Weak form definition
		FEQuadraticV_LinearP fe = new FEQuadraticV_LinearP();
		MathFunc k = C(nu);
		VecMathFunc U = uk;
		MathFunc fc = C0;
		VecMathFunc f = new SpaceVectorFunction(C0, C0);
		VecWeakForm wf = new VecWeakForm(fe,
				(u, v) -> {
					VecMathFunc grad_u1 = grad(u[1],"x","y");
					VecMathFunc grad_v1 = grad(v[1],"x","y");
					VecMathFunc grad_u2 = grad(u[2],"x","y");
					VecMathFunc grad_v2 = grad(v[2],"x","y");
					MathFunc div_u = u[1].diff("x")+u[2].diff("y");
					MathFunc div_v = v[1].diff("x")+v[2].diff("y");
					
					return k * grad_u1.dot(grad_v1)  //   (v1_x,k*u1_x) + (v1_y,k*u1_y)
						 + k * grad_u2.dot(grad_v2)  // + (v2_x,k*u2_x) + (v2_y,k*u2_y)
						 + U.dot(grad_u1)*v[1]      // + (U1*u1_x,v1)+(U2*u1_y,v1)
						 + U.dot(grad_u2)*v[2]      // + (U1*u2_x,v2)+(U2*u2_y,v2)
						 + fc * (u[1]*v[1]+u[2]*v[2])// + c*[(u1,v1)+(u2,v2)]
						 - div_v*u[3]                // - (v1_x+v2_y,p) //where p=u[3]
						 + v[3]*div_u;               // + (q,u1_x+u2_y) //where q=v[3]
				},
				(v)-> v[1]*f[1] + v[2]*f[2]);
		wf.compile();

		VecMathFunc d = new SpaceVectorFunction(2);
		d.set(1, C0);
		d.set(2, C0);
		VecMathFunc normal = new SpaceVectorFunction(2); //normal vector
		normal.set(1, new UserDefFunc() {
			//@Override
			public double apply(AssembleParam ap, double... args) {
				Element be = ap.element;
				Edge edge = (Edge)be.getGeoEntity();
				Vector n = edge.getNormVector();
				return -1.0*n[1];
			}
		});
		normal.set(2, new UserDefFunc() {
			//@Override
			public double apply(AssembleParam ap, double... args) {
				Element be = ap.element;
				Edge edge = (Edge)be.getGeoEntity();
				Vector n = edge.getNormVector();
				return -1.0*n[2];
			}
		});
		VecFiniteElement bfe = fe.getBoundaryFE();
		VecWeakForm wfb = new VecWeakForm(bfe,
				(u, v) -> d[1]*u[1]*v[1] + d[2]*u[2]*v[2],
				(v)    -> v[3]*(normal[1]*v[1]+normal[2]*v[2])
				);
		wfb.compile();
		
		//Define block stiff matrix and block load vector before assembly
		int vvfDim = 3;
		int[] dims = new int[vvfDim];
		for(int vvfIdx=1;vvfIdx<=vvfDim;vvfIdx++) {
			dims[vvfIdx-1] = fe.getNumberOfDOFs(mesh, vvfIdx);
		}
		SparseBlockMatrix stiff = new SparseBlockMatrix(vvfDim,vvfDim);
		SparseBlockVector load = new SparseBlockVector(vvfDim);
		for(int i=1;i<=vvfDim;i++) {
			for(int j=1;j<=vvfDim;j++) {
				stiff.setBlock(i, j, new SparseMatrixRowMajor(dims[i-1],dims[j-1]));
			}
			load.setBlock(i, new SparseVectorHashMap(dims[i-1]));
		}

		BasicVecAssembler assembler = new BasicVecAssembler(mesh, wf);
		//the block matrix and vector are used as normal matrix and vector
		assembler.assembleGlobal(stiff, load); 
		
		// Use BasicAssembler to assemble boundary elements
		BasicVecAssembler boundaryAssembler = new BasicVecAssembler(mesh, wfb);
		for(Element e : mesh.getElementList()) {
			for(Element be : e.getBorderElements()) {
				//Check node type
				////TODO how to check boundary type beside checking node type???
				int nBeDOF = bfe.getNumberOfDOFs();
				for(int i=0; i<nBeDOF; i++) {
					int nVVFIdx = bfe.getVVFComponentIndex(i+1);
					NodeType nodeType = be.getBorderNodeType(nVVFIdx);
					if(nodeType == NodeType.Neumann || nodeType == NodeType.Robin) {
						boundaryAssembler.assembleGlobal(be, stiff, load);
					}
				}
			}
		}
		
		Utils.imposeDirichletCondition(stiff, load, fe, mesh, diri);
		
		Element e = mesh.getElementList().at(1);
		printData(mesh, fe, e, stiff, load);
		
		
		SchurComplementStokesSolver solver = 
			new SchurComplementStokesSolver(stiff,load);
		//solver.setCGInit(0.5);
		//solver.debug = true;
		return solver.solve2D();
	}

	public void run(int testCaseNo) {
		//read mesh and set boundary condiitons
		init(testCaseNo);
		//initial velocity
		SpaceVectorFunction uk = new SpaceVectorFunction(2);
		uk.set(1, C0);
		uk.set(2, C0);
		//solution u
		SparseBlockVector u = null;
		for(int iter=0; iter<this.maxNonlinearIter; iter++) {
			//solve the system with uk as the convection velocity
			u = nonlinearIterSteady(iter, uk);
			//Compute norm of delta_u (not including delta_v)
			int dim = u.getBlock(1).getDim();
			SparseVector delta_u = new SparseVectorHashMap(dim);
			for(int i=1;i<=dim;i++)
				delta_u.set(i, 
						u.getBlock(1).get(i)-uk.get(1).apply(new Variable().setIndex(i)));

//			uk.set(1, new Vector2MathFunc(u.getBlock(1), mesh, "x","y"));
//			uk.set(2, new Vector2MathFunc(u.getBlock(2), mesh, "x","y"));
			SparseVector u1 = u.getBlock(1);
			SparseVector u2 = u.getBlock(2);
			
			uk.set(1, new Vector2MathFunc(u.getBlock(1)));
			uk.set(2, new Vector2MathFunc(u.getBlock(2)));

			System.out.println("Iter="+iter+" Error Norm2 (||u1_k+1 - u1_k||) = "+delta_u.norm2());
			
			Tools.plotVector(mesh, outputFolder, String.format("%s_uv_%02d.dat",file,iter), 
					u.getBlock(1), u.getBlock(2));
			Tools.plotVector(meshOld, outputFolder, String.format("%s_p_%02d.dat",file,iter), 
					Tools.valueOnElement2Node(mesh, u.getBlock(3)));
			if(delta_u.norm2() < this.nonlinearError) {
				break;
			}
		}
	}
	
	public void printData(Mesh mesh, VecFiniteElement fe, Element e, 
			Matrix stiff, Vector load) {
		for(int i=1; i<=fe.getNumberOfDOFs(); i++) {
			int gi = fe.getGlobalIndex(mesh, e, i);
			for(int j=1; j<=fe.getNumberOfDOFs(); j++) {
				int gj = fe.getGlobalIndex(mesh, e, j);
				System.out.print(String.format("%.6f", stiff.get(gi, gj))+"  ");
			}
			System.out.println();
		}
		for(int j=1; j<=fe.getNumberOfDOFs(); j++) {
			int gj = fe.getGlobalIndex(mesh, e, j);
			System.out.println(String.format("%.6f", load.get(gj)));
		}
		
	}
	
	public static void main(String[] args) {
		Ex11_NavierStokesBox NSB = new Ex11_NavierStokesBox();
		System.out.println("mu="+NSB.nu);
		System.out.println("maxNonlinearIter="+NSB.maxNonlinearIter);
		System.out.println("nonlinearError="+NSB.nonlinearError);
		NSB.run(1);
	}
}
