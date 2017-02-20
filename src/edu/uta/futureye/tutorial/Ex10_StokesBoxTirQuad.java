package edu.uta.futureye.tutorial;

import static edu.uta.futureye.function.FMath.C0;
import static edu.uta.futureye.function.FMath.C1;
import static edu.uta.futureye.function.FMath.grad;

import java.util.HashMap;

import edu.uta.futureye.algebra.SparseBlockMatrix;
import edu.uta.futureye.algebra.SparseBlockVector;
import edu.uta.futureye.algebra.SparseMatrixRowMajor;
import edu.uta.futureye.algebra.SparseVectorHashMap;
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
import edu.uta.futureye.core.intf.VecFiniteElement;
import edu.uta.futureye.function.MultiVarFunc;
import edu.uta.futureye.function.UserDefFunc;
import edu.uta.futureye.function.basic.SpaceVectorFunction;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.function.intf.VecMathFunc;
import edu.uta.futureye.io.MeshReader;
import edu.uta.futureye.lib.assembler.AssembleParam;
import edu.uta.futureye.lib.assembler.BasicVecAssembler;
import edu.uta.futureye.lib.element.FEQuadraticV_LinearP;
import edu.uta.futureye.lib.weakform.VecWeakForm;
import edu.uta.futureye.util.Constant;
import edu.uta.futureye.util.Utils;
import edu.uta.futureye.util.container.ElementList;
import edu.uta.futureye.util.container.ObjIndex;

/**
 * Problem:
 *   -\Nabla{k*\Nabla{\vec{u}} + \Nabla{p} = \vec{f}
 *   div{\vec{u}} = 0
 * 
 * Each dim:
 *   -k*(u1_xx+u1_yy) + p_x = f1
 *   -k*(u2_xx+u2_yy) + p_y = f2
 *   u1_x+u2_y              = 0

 * Weak form:
 *   find \vec{u} \in H_0^1(div;\Omega), p \in L_2(\Omega)
 *   such that, for all \vec{v} \in H_0^1(div;\Omega), q \in L_2(\Omega)
 *   
 *   (\Nabla{\vec{v}},k*\Nabla{\vec{u}}) - (div{\vec{v}},p) 
 *                   + (q,div{\vec{u}}) = (\vec{v},\vec{f})
 *
 *   (v1_x,k*u1_x) + (v1_y,k*u1_y) + (v2_x,k*u2_x) + (v2_y,k*u2_y) 
 *                   - (v1_x+v2_y,p) + (q,u1_x+u2_y) = (v1*f1+v2*f2)      
 *
 * where
 *   \vec{u}=(u1,u2): velocity vector field    
 *   \vec{f}=(f1,f2): body force
 *   
 * @author liuyueming
 *
 */
public class Ex10_StokesBoxTirQuad {
	public static String outputFolder = ".";
	
	public static void box() {
		//Read a triangle mesh from an input file
		//[-3,3]*[-3,3]
		String file = "grids/stokes_box";
		
		MeshReader reader = new MeshReader(file+".grd");
		Mesh mesh = reader.read2DMesh();
		mesh.nVertex = mesh.getNodeList().size();
		
		//Add nodes for quadratic element
		for(Element e : mesh.getElementList()) {
			e.adjustVerticeToCounterClockwise();
			int nNode = e.nodes.size();
			for(EdgeLocal edge : e.edges()) {
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
		
		//Geometry relationship
		mesh.computeNodeBelongsToElements();
		
		ElementList eList = mesh.getElementList();
		//NodeList nodes = mesh.getNodeList();
		
		for(int i=1;i<=eList.size();i++) {
			System.out.println(i+"  " + eList.at(i));
		}
		
		//Mark border type
		HashMap<NodeType, MathFunc> mapNTF_uv = new HashMap<NodeType, MathFunc>();
		mapNTF_uv.put(NodeType.Dirichlet, null);
		
		HashMap<NodeType, MathFunc> mapNTF_p = new HashMap<NodeType, MathFunc>();
		mapNTF_p.put(NodeType.Neumann, null);
		
		mesh.markBorderNode(new ObjIndex(1,2),mapNTF_uv);
		mesh.markBorderNode(3,mapNTF_p);

		//Use element library to assign degree of freedom (DOF) to element
		for(int i=1;i<=eList.size();i++) {
			System.out.println(i+"  " + eList.at(i));
		}

		// Weak form definition
		FEQuadraticV_LinearP fe = new FEQuadraticV_LinearP();
		MathFunc k = C1;
		VecMathFunc f = new SpaceVectorFunction(C0, C0);
		VecWeakForm wf = new VecWeakForm(fe,
				(u, v) -> k * grad(u[1],"x","y" ).dot(grad(v[1],"x","y")) //   (v1_x,k*u1_x) + (v1_y,k*u1_y)
						+ k * grad(u[2],"x","y" ).dot(grad(v[2],"x","y")) // + (v2_x,k*u2_x) + (v2_y,k*u2_y) 
						- (v[1].diff("x")+v[2].diff("y"))*u[3]            // - (v1_x+v2_y,p) //where p=u[3]
						+ v[3]*(u[1].diff("x")+u[2].diff("y")),           // + (q,u1_x+u2_y) //where q=v[3]
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

		//Define block stiff matrix and block load vector
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
		assembler.assembleGlobal(stiff, load);
		
		
//		// Use BasicAssembler to assemble boundary elements
//		BasicVecAssembler boundaryAssembler = new BasicVecAssembler(mesh, wfb);
//		for(Element e : eList) {
//			for(Element be : e.getBorderElements()) {
//				//Check node type
//				////TODO how to check boundary type beside checking node type???
//				int nBeDOF = bfe.getNumberOfDOFs();
//				for(int i=0; i<nBeDOF; i++) {
//					int nVVFIdx = bfe.getVVFComponentIndex(i+1);
//					NodeType nodeType = be.getBorderNodeType(nVVFIdx);
//					if(nodeType == NodeType.Neumann || nodeType == NodeType.Robin) {
//						boundaryAssembler.assembleGlobal(be, stiff, load);
//					}
//				}
//			}
//		}
		System.out.println("load norm2="+load.norm2());
		for (int i = 1; i <= 3; i++) {
			for (int j = 1; j <= 3; j++) {
				System.out.println("getNonZeroNumber("+i+","+j+")="+stiff.getBlock(i, j).getNonZeroNumber()+"   ");
			}
		}

		VecMathFunc diri = new SpaceVectorFunction(3);
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
		diri.set(2, C0);
		diri.set(3, C0);

		Utils.imposeDirichletCondition(stiff, load, fe, mesh, diri);
		
		//print some values for debug purpose
		int NB = mesh.getNodeList().size()*2+10;
		int NE = mesh.getNodeList().size()*2+20;
		for (int i = NB; i <= NE; i++) {
			for (int j = NB; j <= NE; j++) {
				System.out.print(stiff.get(i, j)+"     ");
			}
			System.out.println();
		}
		for (int i = NB; i <= NE; i++) {
			System.out.println(load.get(i));
		}
//		
		//TODO how to assemble SparseBlock Matrix and Vector???
		SchurComplementStokesSolver solver = 
			new SchurComplementStokesSolver(stiff,load);
		
		SparseBlockVector u = solver.solve2D();
		
		//没有指定压强，不同的求解器可能会差一个常数
		System.out.println("u=");
		for(int i=1;i<=u.getDim();i++) {
			System.out.println(String.format("%.3f", u.get(i)));
		}
//		Tools.plotVector(mesh, outputFolder, String.format("%s_uv.dat",file), 
//				u.getBlock(1), u.getBlock(2));
	}
	
	public static void main(String[] args) {
		box();
	}
}
