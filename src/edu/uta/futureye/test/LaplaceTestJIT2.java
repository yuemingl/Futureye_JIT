package edu.uta.futureye.test;

import static edu.uta.futureye.function.FMath.*;

import java.util.HashMap;
import java.util.Map;

import edu.uta.futureye.algebra.SparseMatrixRowMajor;
import edu.uta.futureye.algebra.SparseVectorHashMap;
import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.SparseMatrix;
import edu.uta.futureye.algebra.intf.SparseVector;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.algebra.solver.external.SolverJBLAS;
import edu.uta.futureye.bytecode.CompiledFunc;
import edu.uta.futureye.core.DOF;
import edu.uta.futureye.core.DOFOrder;
import edu.uta.futureye.core.Edge;
import edu.uta.futureye.core.EdgeLocal;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Face;
import edu.uta.futureye.core.FaceLocal;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.core.NodeType;
import edu.uta.futureye.core.Volume;
import edu.uta.futureye.core.geometry.GeoEntity;
import edu.uta.futureye.core.geometry.GeoEntity2D;
import edu.uta.futureye.core.geometry.GeoEntity3D;
import edu.uta.futureye.core.intf.WeakForm.ItemType;
import edu.uta.futureye.function.AbstractSimpleMathFunc;
import edu.uta.futureye.function.FMath;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.basic.FX;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.io.MeshReader;
import edu.uta.futureye.io.MeshWriter;
import edu.uta.futureye.lib.assembler.AssemblerScalar;
import edu.uta.futureye.lib.element.FELinearTriangle;
import edu.uta.futureye.lib.shapefun.SFLinearLocal2DRS;
import edu.uta.futureye.lib.weakform.WeakFormLaplace2D;
import edu.uta.futureye.util.FutureyeException;
import edu.uta.futureye.util.container.DOFList;
import edu.uta.futureye.util.container.ElementList;
import edu.uta.futureye.util.container.NodeList;
import edu.uta.futureye.util.container.VertexList;


/**
 * This file expression LHS by defining composite area coordinate variable r(x,y),s(x,y) with symbolic expression
 * <blockquote><pre>
 * Problem:
 *   -\Delta{u} = f
 *   u(x,y)=0, (x,y) \in \partial{\Omega}
 * where
 *   \Omega = [-3,3]*[-3,3]
 *   f = -2*(x^2+y^2)+36
 * Solution:
 *   u = (x^2-9)*(y^2-9)
 * </blockquote></pre>
 * 
 * @author liuyueming
 */
public class LaplaceTestJIT2 {
	public static double eps = 1e-5;
	static double[] triW = {
		0.06296959,
	    0.06619708,
	    0.06296959,
	    0.06619708,
	    0.06296959,
	    0.06619708,
	    0.11250000
	};
	static double[] triR = {
	        0.10128651,
	        0.47014206,
	        0.79742699,
	        0.47014206,
	        0.10128651,
	        0.05971587,
	        0.33333333
	    };
	static double[] triS = {
			0.10128651,
			0.05971587,
			0.10128651,
			0.47014206,
			0.79742699,
			0.47014206,
			0.33333333
		};
	public static double intOnTriangleRefElement(CompiledFunc integrand, double[] params, int paramsStart, int order) {
		double rlt = 0.0;
		if(order == 2) {
			params[paramsStart] = 0.333333333333333;
			params[paramsStart+1] = 0.333333333333333;
			params[paramsStart+2] = 0.333333333333333;
			rlt = 0.5*integrand.apply(params);
		} else if(order == 3) {
			params[paramsStart] = 0.5; params[paramsStart+1] = 0.5; params[paramsStart+2] = 0.0; 
			double pv1 = integrand.apply(params);
			params[paramsStart] = 0.0; params[paramsStart+1] = 0.5; params[paramsStart+2] = 0.5; 
			double pv2 = integrand.apply(params);
			params[paramsStart] = 0.5; params[paramsStart+1] = 0.0; params[paramsStart+2] = 0.5; 
			double pv3 = integrand.apply(params);
			rlt = 0.5*0.333333333333333*(pv1 + pv2 + pv3);
		} else if(order == 4) {
			double w123 = 25.0/48.0;
			double w4 = -27.0/48.0;
			
			params[paramsStart] = 0.6; params[paramsStart+1] = 0.2; params[paramsStart+2] = 0.2; 
			double pv1 = integrand.apply(params);
			params[paramsStart] = 0.2; params[paramsStart+1] = 0.6; params[paramsStart+2] = 0.2; 
			double pv2 = integrand.apply(params);
			params[paramsStart] = 0.2; params[paramsStart+1] = 0.2; params[paramsStart+2] = 0.6; 
			double pv3 = integrand.apply(params);
			params[paramsStart] = 0.333333333333333; params[paramsStart+1] = 0.333333333333333; params[paramsStart+2] = 0.333333333333333; 
			double pv4 = 0.5*integrand.apply(params);
			
			rlt = 0.5*w123*(pv1 + pv2 + pv3) + w4*pv4;
		} else if(order == 5) {
			for(int i=0;i<7;i++) {
				params[paramsStart]   = triR[i]; 
				params[paramsStart+1] = triS[i]; 
				params[paramsStart+2] = 1.0-triR[i]-triS[i];
				rlt += triW[i]*integrand.apply(params);
			}
		}
		return rlt;
	}
	
	public Mesh mesh;
	public Vector u;
	
	public static class TriAreaCoordR extends AbstractSimpleMathFunc {
		MathFunc jac;
		FX x1 = new FX("x1");
		FX x2 = new FX("x2");
		FX x3 = new FX("x3");
		FX y1 = new FX("y1");
		FX y2 = new FX("y2");
		FX y3 = new FX("y3");
		public TriAreaCoordR() {
			super("r", "r");
		}
		public void setJac(MathFunc jac) {
			this.jac = jac;
		}

		@Override
		public double apply(double... args) {
			return 0;
		}
		
//		r_x = (y2-y3)/jac;
//		r_y = (x3-x2)/jac;
		@Override
		public MathFunc diff(String varName) {
			if(varName.equals("r"))
				return FMath.C1;
			if(varName.equals("x"))
				return (y2-y3)/jac;
			else if(varName.equals("y"))
				return (x3-x2)/jac;
			else
				return FMath.C0;
		}
		public String toString() {
			return this.fName;
		}
		public String getExpr() {
			return this.varName;
		}
	}
	public static class TriAreaCoordS extends AbstractSimpleMathFunc {
		MathFunc jac;
		FX x1 = new FX("x1");
		FX x2 = new FX("x2");
		FX x3 = new FX("x3");
		FX y1 = new FX("y1");
		FX y2 = new FX("y2");
		FX y3 = new FX("y3");
		public TriAreaCoordS() {
			super("s", "s");
		}
		public void setJac(MathFunc jac) {
			this.jac = jac;
		}

		@Override
		public double apply(double... args) {
			return 0;
		}
		
//		s_x = (y3-y1)/jac;
//		s_y = (x1-x3)/jac;
		@Override
		public MathFunc diff(String varName) {
			if(varName.equals("s"))
				return FMath.C1;
			if(varName.equals("x"))
				return (y3-y1)/jac;
			else if(varName.equals("y"))
				return (x1-x3)/jac;
			else
				return FMath.C0;
		}
		
		public String toString() {
			return this.fName;
		}
		public String getExpr() {
			return this.varName;
		}
	}
	public void run() {
        //1.Read in a triangle mesh from an input file with
        //  format ASCII UCD generated by Gridgen
        MeshReader reader = new MeshReader("triangle.grd");
        Mesh mesh = reader.read2DMesh();
        //Compute geometry relationship between nodes and elements
        mesh.computeNodeBelongsToElements();

        //2.Mark border types
        HashMap<NodeType, MathFunc> mapNTF =
                new HashMap<NodeType, MathFunc>();
        mapNTF.put(NodeType.Dirichlet, null);
        mesh.markBorderNode(mapNTF);

        //3.Use element library to assign degrees of
        //  freedom (DOF) to element
        ElementList eList = mesh.getElementList();
        FELinearTriangle feLT = new FELinearTriangle();
        for(int i=1;i<=eList.size();i++)
            feLT.assignTo(eList.at(i));

		//Construct a function with the coordinate of points in an element as parameters
		String[] argsOrder = new String[]{"x1","x2","x3","y1","y2","y3","r","s","t"};
		FX x1 = new FX("x1");
		FX x2 = new FX("x2");
		FX x3 = new FX("x3");
		FX y1 = new FX("y1");
		FX y2 = new FX("y2");
		FX y3 = new FX("y3");
		MathFunc fx = x1*r + x2*s + x3*(1-r-s);
		MathFunc fy = y1*r + y2*s + y3*(1-r-s);
		Map<String, MathFunc> map = new HashMap<String, MathFunc>();
		map.put("x", fx);
		map.put("y", fy);
		
		//             (r[0] r[1])   (x_r, x_s)
		// 2D JacMat = (r[2] r[3]) = (y_r, y_s)
		//jac changes with element, define the expression for jac with linear element
		MathFunc jac = fx.diff("r")*fy.diff("s") - fy.diff("r")*fx.diff("s");

        //Right hand side(RHS): f = -2*(x^2+y^2)+36
        MathFunc f = -2*(x*x+y*y)+36;
		MathFunc ff = f.compose(map);
		
		//4.Weak form
		Element ee = eList.at(1); //a template element
		DOFList DOFs = ee.getAllDOFList(DOFOrder.NEFV);
		int nDOFs = DOFs.size();
		
		//形函数计算需要和单元关联
		for(int i=1;i<=nDOFs;i++) {
			DOFs.at(i).getSSF().assignElement(ee);
		}
		MathFunc[][] lhs = new MathFunc[nDOFs][nDOFs];
		MathFunc[] rhs = new MathFunc[nDOFs];
		for(int j=0; j<nDOFs; j++) {
			DOF dofI = DOFs.at(j+1);
			MathFunc v = dofI.getSSF();
//			for(int i=0; i<nDOFs; i++) {
//				DOF dofJ = DOFs.at(i+1);
//				MathFunc u = dofJ.getSSF();
//				//lhs[j][i] = (u.diff("x")*v.diff("x")+u.diff("y")*v.diff("y"))*jac;
//				//lhs[j][i] = (grad(u,"x","y").dot(grad(v,"x","y")))*jac;
//				lhs[j][i].setName("LHS"+i+""+j);
//			}
			rhs[j] = v*ff*jac;//80.966
			rhs[j].setName("RHS"+j);
		}
		
//		r_x = (y2-y3)/jac;
//		r_y = (x3-x2)/jac;
//		s_x = (y3-y1)/jac;
//		s_y = (x1-x3)/jac;
//		double[][] lhs = { //SymJava Example6
//				{rx*rx + ry*ry, rx*sx + ry*sy, rx*tx+ry*ty},
//				{sx*rx + sy*ry, sx*sx + sy*sy, sx*tx+sy*ty},
//				{tx*rx + ty*ry, tx*sx + ty*sy, tx*tx+ty*ty},
//		};
		
//		lhs[0][0] = ((y2-y3)/jac*(y2-y3)/jac             + (x3-x2)/jac*(x3-x2)/jac)*jac; --simplify-->
		//singleton?
		TriAreaCoordR rr = new TriAreaCoordR();
		rr.setJac(jac);
		TriAreaCoordS ss = new TriAreaCoordS();
		ss.setJac(jac);
		MathFunc[] sf = new MathFunc[3];
		sf[0] = rr;
		sf[1] = ss;
		sf[2] = 1-rr-ss;
		for(int j=0; j<nDOFs; j++) {
			MathFunc v = sf[j];
			System.out.println(">>>"+sf[j]);
			for(int i=0; i<nDOFs; i++) {
				MathFunc u = sf[i];
				//lhs[j][i] = (u.diff("x")*v.diff("x")+u.diff("y")*v.diff("y"))*jac;
				lhs[j][i] = (grad(u,"x","y").dot(grad(v,"x","y")))*jac;
				lhs[j][i].setName("LHS"+i+""+j);
			}
			rhs[j] = v*ff*jac;//81.580 
			rhs[j].setName("RHS"+j);
		}
		
//		lhs[0][0] = ((y2-y3)*(y2-y3)             + (x3-x2)*(x3-x2))/jac;
//		lhs[0][1] = ((y2-y3)*(y3-y1)             + (x3-x2)*(x1-x3))/jac;
//		lhs[0][2] = ((y2-y3)*(-(y2-y3)-(y3-y1))  + (x3-x2)*(-(x3-x2)-(x1-x3)))/jac;
//
//		lhs[1][0] = ((y3-y1)*(y2-y3)             + (x1-x3)*(x3-x2))/jac;
//		lhs[1][1] = ((y3-y1)*(y3-y1)             + (x1-x3)*(x1-x3))/jac;
//		lhs[1][2] = ((y3-y1)*(-(y2-y3)-(y3-y1))  + (x1-x3)*(-(x3-x2)-(x1-x3)))/jac;
//		
//		lhs[2][0] = ((-(y2-y3)-(y3-y1))*(y2-y3)             + (-(x3-x2)-(x1-x3))*(x3-x2))/jac;
//		lhs[2][1] = ((-(y2-y3)-(y3-y1))*(y3-y1)             + (-(x3-x2)-(x1-x3))*(x1-x3))/jac;
//		lhs[2][2] = ((-(y2-y3)-(y3-y1))*(-(y2-y3)-(y3-y1))  + (-(x3-x2)-(x1-x3))*(-(x3-x2)-(x1-x3)))/jac;
//		for(int j=0; j<nDOFs; j++) {
//			for(int i=0; i<nDOFs; i++) {
//				lhs[j][i].setName("LHS"+i+""+j);
//			}
//			rhs[j].setName("RHS"+j);
//		}


		CompiledFunc[][] clhs = new CompiledFunc[nDOFs][nDOFs];
		CompiledFunc[] crhs = new CompiledFunc[nDOFs];
		for(int j=0; j<nDOFs; j++) {
			for(int i=0; i<nDOFs; i++) {
				clhs[j][i] = lhs[j][i].compile(argsOrder);
			}
			crhs[j] = rhs[j].compile(argsOrder);
		}
		
		//5.Assembly process
		double[][] A = new double[nDOFs][nDOFs];
		double[] b = new double[nDOFs];
		double[] params = new double[argsOrder.length];
		int dim = mesh.getNodeList().size();
		SparseMatrix stiff = new SparseMatrixRowMajor(dim,dim);;
		SparseVector load = new SparseVectorHashMap(dim);
		
		long start = System.currentTimeMillis();
		for(Element e : eList) {
			//e.adjustVerticeToCounterClockwise();

			DOFs = e.getAllDOFList(DOFOrder.NEFV);
			double[] coords = e.getNodeCoords();
			System.arraycopy(coords, 0, params, 0, coords.length);

			for(int j=0; j<nDOFs; j++) {
				for(int i=0; i<nDOFs; i++) {
					A[j][i] = intOnTriangleRefElement(clhs[j][i], params, coords.length, 3);
				}
				b[j] = intOnTriangleRefElement(crhs[j], params, coords.length, 3);
			}
			
			for(int j=0;j<nDOFs;j++) {
				DOF dofI = DOFs.at(j+1);
				int nGlobalRow = dofI.getGlobalIndex();
				for(int i=0;i<nDOFs;i++) {
					DOF dofJ = DOFs.at(i+1);
					int nGlobalCol = dofJ.getGlobalIndex();
					stiff.add(nGlobalRow, nGlobalCol, A[j][i]);
				}
				//Local load vector
				load.add(nGlobalRow, b[j]);
			}
		}
		System.out.println("Aassembly time: "+(System.currentTimeMillis()-start)+"ms");

		//Boundary condition
		imposeDirichletCondition(stiff, load, mesh, C0);
		
/*        for(int i=1; i<=stiff.getRowDim(); i++) {
        	for(int j=1; j<=stiff.getRowDim(); j++) {
        		System.out.print(stiff.get(i,j)+" ");
        	}
        	System.out.println();
        }
        for(int i=1; i<=load.getDim(); i++) {
        	System.out.print(load.get(i)+" ");
        } */

        //6.Solve linear system
        SolverJBLAS solver = new SolverJBLAS();
        Vector u = solver.solveDGESV(stiff, load);
        System.out.println("u=");
        for(int i=1;i<=u.getDim();i++)
            System.out.println(String.format("%.3f ", u.get(i)));

        //7.Output results to an Techplot format file
        MeshWriter writer = new MeshWriter(mesh);
        writer.writeTechplot("./tutorial/Laplace2D.dat", u);

        this.mesh = mesh;
        this.u = u;
	}
	
	public static void setDirichlet(Matrix stiff, Vector load, int matIndex, double value) {
		int row = matIndex;
		int col = matIndex;
		stiff.set(row, col, 1.0);
		load.set(row,value);
		for(int r=1;r<=stiff.getRowDim();r++) {
			if(r != row) {
				load.add(r,-stiff.get(r, col)*value);
				stiff.set(r, col, 0.0);
			}
		}
		for(int c=1;c<=stiff.getColDim();c++) {
			if(c != col) {
				stiff.set(row, c, 0.0);
			}
		}
	}
	
	public static void imposeDirichletCondition(Matrix stiff, Vector load, Mesh mesh, MathFunc diri) {
		ElementList eList = mesh.getElementList();
		for(int i=1;i<=eList.size();i++) {
			Element e = eList.at(i);
			DOFList DOFs = e.getAllDOFList(DOFOrder.NEFV);
			for(int j=1;j<=DOFs.size();j++) {
				DOF dof = DOFs.at(j);
				GeoEntity ge = dof.getOwner();
				if(ge instanceof Node) {
					Node n = (Node)ge;
					if(n.getNodeType() == NodeType.Dirichlet) {
						Variable v = Variable.createFrom(diri, n, n.globalIndex); //bugfix 11/27/2013 Variable.createFrom(diri, n, 0);
						setDirichlet(stiff, load, dof.getGlobalIndex(),diri.apply(v));
					}
				} else if(ge instanceof EdgeLocal) {
					//2D单元（面）其中的局部边上的自由度
					EdgeLocal edge = (EdgeLocal)ge;
					if(edge.getBorderType() == NodeType.Dirichlet) {
						//TODO 以边的那个顶点取值？中点？
						//Variable v = Variable.createFrom(fdiri, ?, 0);
					}
					
				} else if(ge instanceof FaceLocal) {
					//3D单元（体）其中的局部面上的自由度
					FaceLocal face = (FaceLocal)ge;
					if(face.getBorderType() == NodeType.Dirichlet) {
						//TODO
					}
				} else if(ge instanceof Edge) {
					//1D单元（线段）上的自由度，其Dirichlet边界用结点来计算推出，而不需要专门标记单元
					VertexList vs = ((GeoEntity2D) ge).getVertices();
					for(int k=1;k<=vs.size();k++) {
						Node n = vs.at(k).globalNode();
						if(NodeType.Dirichlet == n.getNodeType()) {
							Variable v = Variable.createFrom(diri, n, 0);
							setDirichlet(stiff, load, dof.getGlobalIndex(),diri.apply(v));
						}
					}
				} else if(ge instanceof Face) {
					//2D单元（面）上的自由度，其Dirichlet边界用结点来计算推出，而不需要专门标记单元
					
					VertexList vs = ((GeoEntity2D) ge).getVertices();
					for(int k=1;k<=vs.size();k++) {
						Node n = vs.at(k).globalNode();
						if(NodeType.Dirichlet == n.getNodeType()) {
							Variable v = Variable.createFrom(diri, n, 0);
							setDirichlet(stiff, load, dof.getGlobalIndex(),diri.apply(v));
						}
					}
				} else if(ge instanceof Volume) {
					//3D单元（体）上的自由度，其Dirichlet边界用结点来计算推出，而不需要专门标记单元
					VertexList vs = ((GeoEntity3D) ge).getVertices();
					for(int k=1;k<=vs.size();k++) {
						Node n = vs.at(k).globalNode();
						if(NodeType.Dirichlet == n.getNodeType()) {
							Variable v = Variable.createFrom(diri, n, 0);
							setDirichlet(stiff, load, dof.getGlobalIndex(),diri.apply(v));
						}
					}
				}
			}
		}
	}	
    public static void main(String[] args) {
    	LaplaceTestJIT2 ex1 = new LaplaceTestJIT2();
    	ex1.run();
    }
}
