package edu.uta.futureye.test;

import static edu.uta.futureye.function.FMath.*;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
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
import edu.uta.futureye.function.FMath;
import edu.uta.futureye.function.SingleVarFunc;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.basic.FX;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.function.operator.FOIntegrate;
import edu.uta.futureye.io.MeshReader;
import edu.uta.futureye.io.MeshWriter;
import edu.uta.futureye.lib.assembler.AssemblerScalar;
import edu.uta.futureye.lib.element.FELinearTriangle;
import edu.uta.futureye.lib.shapefun.SFLinearLocal2DRS;
import edu.uta.futureye.lib.weakform.WeakFormLaplace2D;
import edu.uta.futureye.util.FutureyeException;
import edu.uta.futureye.util.Utils;
import edu.uta.futureye.util.container.DOFList;
import edu.uta.futureye.util.container.ElementList;
import edu.uta.futureye.util.container.NodeList;
import edu.uta.futureye.util.container.VertexList;


/**
 * weakform and lambda expression are used
 * for This file LHS anf RHS

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
public class LaplaceTestJIT3 {
	
	public Mesh mesh;
	public Vector u;
	
	public static class TriAreaCoordR extends SingleVarFunc {
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
	public static class TriAreaCoordS extends SingleVarFunc {
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
	
	public interface LHSExpr {
		MathFunc apply(MathFunc u, MathFunc v);
	}
	
	public interface RHSExpr {
		MathFunc apply(MathFunc v);
	}
	
	public static class FELinearTriangleT {
		//Construct a function with the coordinate of points in an element as parameters
		String[] argsOrder = new String[]{"x1","x2","x3","y1","y2","y3","r","s","t"};
		FX x1 = new FX("x1");
		FX x2 = new FX("x2");
		FX x3 = new FX("x3");
		FX y1 = new FX("y1");
		FX y2 = new FX("y2");
		FX y3 = new FX("y3");
		MathFunc fx;
		MathFunc fy;
		Map<String, MathFunc> map;
		MathFunc jac;
		MathFunc[][] matLHS;
		MathFunc[] vecRHS;
		public int nDOFs = 3;
		
		public FELinearTriangleT() {
			fx = x1*r + x2*s + x3*(1-r-s);
			fy = y1*r + y2*s + y3*(1-r-s);
			map = new HashMap<String, MathFunc>();
			map.put("x", fx);
			map.put("y", fy);
			//             (r[0] r[1])   (x_r, x_s)
			// 2D JacMat = (r[2] r[3]) = (y_r, y_s)
			//jac changes with element, define the expression for jac with linear element
			jac = fx.diff("r")*fy.diff("s") - fy.diff("r")*fx.diff("s");
			matLHS = new MathFunc[nDOFs][nDOFs];
			vecRHS = new MathFunc[nDOFs];
		}
		
		public void makeWeakForm(LHSExpr lhsExpr, RHSExpr rhsExpr) {
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
					matLHS[j][i] = lhsExpr.apply(u, v).compose(map)*jac;
					matLHS[j][i].setName("LHS"+i+""+j);
				}
				vecRHS[j] = rhsExpr.apply(v).compose(map)*jac;
				vecRHS[j].setName("RHS"+j);
			}
		}
		
		CompiledFunc[][] clhs = new CompiledFunc[nDOFs][nDOFs];
		CompiledFunc[] crhs = new CompiledFunc[nDOFs];
		
		public void compileWeakForm() {
			clhs = new CompiledFunc[nDOFs][nDOFs];
			crhs = new CompiledFunc[nDOFs];
			for(int j=0; j<nDOFs; j++) {
				for(int i=0; i<nDOFs; i++) {
					clhs[j][i] = matLHS[j][i].compile(argsOrder);
				}
				crhs[j] = vecRHS[j].compile(argsOrder);
			}
		}
		
		public CompiledFunc[][] getCompiledLHS() {
			return clhs;
		}
		
		public CompiledFunc[] getCompiledRHS() {
			return crhs;
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
		FELinearTriangleT fet = new FELinearTriangleT();
		
        //Right hand side(RHS): f = -2*(x^2+y^2)+36
        MathFunc f = -2*(x*x+y*y)+36;

		//4.Weak form
        //The direct implementation of the weak form requires a user write the expression of the weak form inside loops. 
        //Users must responsible to the correctness of loops. 
        //In order to avoid writing the loops by the users, the idea of a template expression of the weak form can be adopted. 
        //By providing by users just the template expression of the weak form, the trial and test functions in the expression 
        //can be replaced to concrete shape functions in the loops of the library code.
        //There are several ways in the implementation of the replacement of the trial and test functions.
        //(1) A straight forward way is using the replace method for symbolic expression to replace the symbol of trial and 
        //test functions by the symbol of shape functions. This way has two drawbacks. First, it is slow since the replacement operation 
        //is actually a string match. Second, symbols could be replace by wrong ones during the string matching without resulting any error messages.
        //(2) Let the users define a weak form function with trial and test functions as parameters. This way provide fast speed and better error checking. 
        //However, the user interface is cumbersome since the users have to define functions for the weak forms.
        //(3) Follow the idea from (2), but using the new feature lambda expression provided by Java 8, the expression of the weak form can be define concisely 
        //with all the advantages of method (2). Specifically, we define two functional interfaces of the weak form builder to accept 
        //the left hand side and right hand side of a weak form by providing two lambda expression by the users.
		//Java 8:
        fet.makeWeakForm(
				(u,v) -> grad(u,"x","y").dot(grad(v,"x","y")), 
				v -> f*v
		);
		fet.compileWeakForm();

		CompiledFunc[][] clhs = fet.getCompiledLHS();
		CompiledFunc[] crhs = fet.getCompiledRHS();
		int nDOFs = fet.nDOFs;
		
		//5.Assembly process
		double[][] A = new double[nDOFs][nDOFs];
		double[] b = new double[nDOFs];
		double[] params = new double[argsOrder.length];
		int dim = mesh.getNodeList().size();
		SparseMatrix stiff = new SparseMatrixRowMajor(dim,dim);
		SparseVector load = new SparseVectorHashMap(dim);
		
		long start = System.currentTimeMillis();
		for(Element e : eList) {
			//e.adjustVerticeToCounterClockwise();

			DOFList DOFs = e.getAllDOFList(DOFOrder.NEFV);
			double[] coords = e.getNodeCoords();
			System.arraycopy(coords, 0, params, 0, coords.length);

			for(int j=0; j<nDOFs; j++) {
				for(int i=0; i<nDOFs; i++) {
					A[j][i] = FOIntegrate.intOnTriangleRefElement(clhs[j][i], params, coords.length, 3);
				}
				b[j] = FOIntegrate.intOnTriangleRefElement(crhs[j], params, coords.length, 3);
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
		Utils.imposeDirichletCondition(stiff, load, mesh, C0);
		
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
	
    public static void main(String[] args) {
    	LaplaceTestJIT3 ex1 = new LaplaceTestJIT3();
    	ex1.run();
    }
}
