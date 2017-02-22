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
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.NodeType;
import edu.uta.futureye.function.MultiVarFunc;
import edu.uta.futureye.function.UserDefFunc;
import edu.uta.futureye.function.basic.SpaceVectorFunction;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.function.intf.VecMathFunc;
import edu.uta.futureye.io.MeshReader;
import edu.uta.futureye.lib.assembler.AssembleParam;
import edu.uta.futureye.lib.assembler.BasicVecAssembler;
import edu.uta.futureye.lib.element.FEBilinearV_ConstantP;
import edu.uta.futureye.lib.weakform.VecWeakForm;
import edu.uta.futureye.util.Constant;
import edu.uta.futureye.util.Utils;
import edu.uta.futureye.util.container.ElementList;
import edu.uta.futureye.util.container.ObjIndex;

/**
 * Stokes problem is defined as
 *  -\nabla\cdot(k\nabla{\vec{u}}) + \nabla{p} = \vec{f}\\
 *  div~{\vec{u}} = 0
 *
 *  or written explicitly in each dimension:
 *   -k*(u1_xx+u1_yy) + p_x = f1
 *   -k*(u2_xx+u2_yy) + p_y = f2
 *   u1_x+u2_y              = 0
 *
 * The weak form is
 *   find \vec{u} \in H_0^1(div;\Omega), p \in L_2(\Omega)
 *   such that, for all \vec{v} \in H_0^1(div;\Omega), q \in L_2(\Omega)
 *   
 *   (\nabla{\vec{v}},k\nabla{\vec{u}}) - (div~{\vec{v}},p) 
 *                   + (q,div~{\vec{u}}) = (\vec{v},\vec{f})
 *   or written explicitly:
 *   (v1_x,k*u1_x) + (v1_y,k*u1_y) + (v2_x,k*u2_x) + (v2_y,k*u2_y) 
 *                 - (v1_x+v2_y,p) + (q,u1_x+u2_y) = (v1*f1+v2*f2)
 * where
 *   \vec{u}=(u1,u2): velocity vector field
 *   \vec{f}=(f1,f2): body force
 *
 */
public class Ex10_StokesBox {
	public static String outputFolder = ".";
	
	public static void box() {
		//Read a rectanlge mesh from an input file
		//[-3,3]*[-3,3]
		String file = "grids/rectangle";
		
		MeshReader reader = new MeshReader(file+".grd");
		Mesh mesh = reader.read2DMesh();
		mesh.nVertex = mesh.getNodeList().size();
		
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

//		//Use element library to assign degree of freedom (DOF) to element
//		for(int i=1;i<=eList.size();i++) {
//			System.out.println(i+"  " + eList.at(i));
//		}

		// Weak form definition
		FEBilinearV_ConstantP fe = new FEBilinearV_ConstantP();
		MathFunc k = C1;
		VecMathFunc f = new SpaceVectorFunction(C0, C0);
		VecWeakForm wf = new VecWeakForm(fe,
				(u, v) ->
					//  (v1_x,k*u1_x) + (v1_y,k*u1_y) 
					//+ (v2_x,k*u2_x) + (v2_y,k*u2_y) 
					//- (v1_x+v2_y,p) 
					//+ (q,u1_x+u2_y)
					//where p=u[3]; q=v[3]
					  k * grad(u[1],"x","y" ).dot(grad(v[1],"x","y"))
					+ k * grad(u[2],"x","y" ).dot(grad(v[2],"x","y"))
					- (v[1].diff("x")+v[2].diff("y"))*u[3]
					+ v[3]*(u[1].diff("x")+u[2].diff("y")),
				(v)-> v[1]*f[1] + v[2]*f[2]
				);
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
		VecWeakForm wfb = new VecWeakForm(fe.getBoundaryFE(),
				(u, v) -> d[1]*u[1]*v[1] + d[2]*u[2]*v[2],
				(v)    -> v[3]*(normal[1]*v[1]+normal[2]*v[2])
				);
		wfb.compile();

		int vvfDim = 3;
		int[] dims = new int[vvfDim];
		for(int vvfIdx=1;vvfIdx<=vvfDim;vvfIdx++) {
			dims[vvfIdx-1] = fe.getNumberOfDOFs(mesh, vvfIdx);
		}
		SparseBlockMatrix stiff = new SparseBlockMatrix(vvfDim,vvfDim);
		SparseBlockVector load = new SparseBlockVector(vvfDim);
		for(int i=1;i<=vvfDim;i++) {
			for(int j=1;j<=vvfDim;j++) {
				stiff.setBlock(i, j, 
						new SparseMatrixRowMajor(dims[i-1],dims[j-1]));
			}
			load.setBlock(i, 
					new SparseVectorHashMap(dims[i-1]));
		}

		BasicVecAssembler assembler = new BasicVecAssembler(mesh, wf);
		assembler.assembleGlobal(stiff, load);
		
		// Use BasicAssembler to assemble boundary elements
		BasicVecAssembler boundaryAssembler = new BasicVecAssembler(mesh, wfb);
		for(Element e : eList) {
			for(Element be : e.getBorderElements()) {
				//Check node type
				////TODO how to check boundary type beside checking node type???
				NodeType nodeType = be.getBorderNodeType();
				if(nodeType == NodeType.Neumann || nodeType == NodeType.Robin) {
					
					boundaryAssembler.assembleGlobal(be, 
							assembler.getGlobalStiffMatrix(), assembler.getGlobalLoadVector());
				}
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
//		load.getBlock(1).print();
//		load.getBlock(2).print();
//		load.getBlock(3).print();
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
	    
		Tools.plotVector(mesh, outputFolder, String.format("%s_uv.dat",file), 
				u.getBlock(1), u.getBlock(2));
	}
	
	public static void main(String[] args) {
		box();
	}
}
