package edu.uta.futureye.tutorial;

import static edu.uta.futureye.function.FMath.*;

import java.util.HashMap;

import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.algebra.solver.external.SolverJBLAS;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.NodeType;
import edu.uta.futureye.core.intf.WeakForm;
import edu.uta.futureye.function.AbstractMathFunc;
import edu.uta.futureye.function.FMath;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.function.intf.ScalarShapeFunction;
import edu.uta.futureye.io.MeshReader;
import edu.uta.futureye.io.MeshWriter;
import edu.uta.futureye.lib.assembler.AssemblerScalar;
import edu.uta.futureye.lib.element.FEBilinearRectangleRegular;
import edu.uta.futureye.lib.element.FELinearTriangle;
import edu.uta.futureye.lib.weakform.WeakFormBuilder;
import edu.uta.futureye.util.container.ElementList;
import edu.uta.futureye.util.container.ObjList;

/**
 * Demo for how to use WeakFormBuilder
 * 
 * @author liuyueming
 *
 */
public class T04UseWeakFormBuilder {
	
	/**
	 * <blockquote><pre>
	 * Solve
	 *   -\Delta{u} = f  in \Omega
	 *   u(x,y) = 0,     on \partial{\Omega}
	 * where
	 *   \Omega = [-3,3]*[-3,3]
	 *   f = -2*(x^2+y^2)+36
	 * Solution:
	 *   u = (x^2-9)*(y^2-9)
	 * </blockquote></pre>   
	 */
	public static void triangleMesh() {
        //1.Read in a triangle mesh from an input file with
        //  format ASCII UCD generated by Gridgen
        MeshReader reader = new MeshReader("triangle.grd");
        Mesh mesh = reader.read2DMesh();
        //Compute geometry relationship of nodes and elements
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

        //4.Weak form. We use WeakFormBuilder to define weak form
		WeakFormBuilder wfb = new WeakFormBuilder() {
			/**
			 * Override this function to define weak form
			 */
			@Override
			public MathFunc makeExpression(Element e, Type type) {
				ScalarShapeFunction u = getScalarTrial();
				ScalarShapeFunction v = getScalarTest();
				MathFunc ff = getParam("f",e);
				switch(type) {
					case LHS_Domain:
						// A(u,v) = u_x*v_x + u_y*v_y
						return grad(u, "x","y").dot(grad(v, "x","y"));
						//return u._d("x").M(v._d("x")) .A (u._d("y").M(v._d("y")));
					case RHS_Domain:
						return ff.M(v);
					default:
						return null;
				}
			}
		};
        //Right hand side(RHS): f = -2*(x^2+y^2)+36
        wfb.addParam("f",x.M(x).A(y.M(y)).M(-2.0).A(36.0));
		
        //5.Assembly process
        AssemblerScalar assembler =
                new AssemblerScalar(mesh, wfb.getScalarWeakForm());
        assembler.assemble();
        Matrix stiff = assembler.getStiffnessMatrix();
        Vector load = assembler.getLoadVector();
        //Boundary condition
        assembler.imposeDirichletCondition(C0);

        //6.Solve linear system
        SolverJBLAS solver = new SolverJBLAS();
        Vector u = solver.solveDGESV(stiff, load);
        System.out.println("u=");
        for(int i=1;i<=u.getDim();i++)
            System.out.println(String.format("%.3f", u.get(i)));

        //7.Output results to an Techplot format file
        MeshWriter writer = new MeshWriter(mesh);
        writer.writeTechplot("./tutorial/UseWeakFormBuilder1.dat", u);
	}
	
	
	/**
	 * <blockquote><pre>
	 * Solve
	 *   -k*\Delta{u} = f  in \Omega
	 *   u(x,y) = 0,       on boundary x=3.0 of \Omega
	 *   u_n + u = 0.01,   on other boundary of \Omega
	 * where
	 *   \Omega = [-3,3]*[-3,3]
	 *   k = 2
	 *   f = -4*(x^2+y^2)+72
	 *   u_n = \frac{\partial{u}}{\partial{n}}
	 *   n: outer unit normal of \Omega
	 * </blockquote></pre>
	 */
	public static void rectangleTest() {
        //1.Read in a rectangle mesh from an input file with
        //  format ASCII UCD generated by Gridgen
		MeshReader reader = new MeshReader("rectangle.grd");
		Mesh mesh = reader.read2DMesh();
		mesh.computeNodeBelongsToElements();
		
        //2.Mark border types
		HashMap<NodeType, MathFunc> mapNTF = new HashMap<NodeType, MathFunc>();
		//Robin type on boundary x=3.0 of \Omega
		mapNTF.put(NodeType.Robin, new AbstractMathFunc("x","y"){
			@Override
			public double apply(Variable v) {
				if(3.0-v.get("x") < 0.01)
					return 1.0; //this is Robin condition
				else
					return -1.0;
			}
		});
		//Dirichlet type on other boundary of \Omega
		mapNTF.put(NodeType.Dirichlet, null);		
		mesh.markBorderNode(mapNTF);

        //3.Use element library to assign degrees of
        //  freedom (DOF) to element
        ElementList eList = mesh.getElementList();
		//FEBilinearRectangle bilinearRectangle = new FEBilinearRectangle();
        //If the boundary of element parallel with coordinate use this one instead.
        //It will be faster than the old one.
		FEBilinearRectangleRegular bilinearRectangle = new FEBilinearRectangleRegular();
        for(int i=1;i<=eList.size();i++)
        	bilinearRectangle.assignTo(eList.at(i));		
		
		
        //4.Weak form. We use WeakFormBuilder to define weak form
		WeakFormBuilder wfb = new WeakFormBuilder() {
			/**
			 * Override this function to define weak form
			 */
			@Override
			public MathFunc makeExpression(Element e, Type type) {
				ScalarShapeFunction u = getScalarTrial();
				ScalarShapeFunction v = getScalarTest();
				//Call param() to get parameters, do NOT define functions here 
				//except for constant functions (or class FC). Because functions 
				//will be transformed to local coordinate system by param()
				MathFunc fk = getParam("k",e);
				MathFunc ff = getParam("f",e);
				switch(type) {
					case LHS_Domain:
						// k*(u_x*v_x + u_y*v_y) in \Omega
						return fk.M(grad(u, "x","y").dot(grad(v, "x","y")));
						//return fk.M( u._d("x").M(v._d("x")) .A (u._d("y").M(v._d("y"))) );
					case LHS_Border:
						// k*u*v on Robin boundary
						return fk.M(u.M(v));
					case RHS_Domain:
						return ff.M(v);
					case RHS_Border:
						return v.M(0.01);
					default:
						return null;
				}
			}
		};
        wfb.addParam("k", C(2.0));
        //Right hand side(RHS): f = -4*(x^2+y^2)+72
        wfb.addParam("f",x.M(x).A(y.M(y)).M(-4.0).A(72.0));
		
        //5.Assembly process
        AssemblerScalar assembler =
                new AssemblerScalar(mesh, wfb.getScalarWeakForm());
        System.out.println("Begin Assemble...");
        assembler.assemble();
        Matrix stiff = assembler.getStiffnessMatrix();
        Vector load = assembler.getLoadVector();
        //Boundary condition
        assembler.imposeDirichletCondition(C0);
        System.out.println("Assemble done!");

        //6.Solve linear system
        SolverJBLAS solver = new SolverJBLAS();
        Vector u = solver.solveDGESV(stiff, load);
        System.out.println("u=");
        for(int i=1;i<=u.getDim();i++)
            System.out.println(String.format("%.3f", u.get(i)));

        //7.Output results to an Techplot format file
        MeshWriter writer = new MeshWriter(mesh);
        writer.writeTechplot("./tutorial/UseWeakFormBuilder2.dat", u);
	}

	
	/**
	 * <blockquote><pre>
	 * Solve
	 *   -k*\Delta{u} + c*u = f  in \Omega
	 *   u(x,y) = u0,     on \Gamma_D
	 *   d*u + k*u_n = g  on \Gamma_N
	 *   
	 * where
	 *   \Omega = [-3,3]*[-3,3]
	 *   \Gamma = \partial{\Omega} (Boundary of \Omega)
	 *   f = -2*(x^2+y^2)+36
	 *   k = 0.1
	 *   c = 1.0
	 *   d = 1.0
	 *   g = x^2+y^2
	 *   
	 * </blockquote></pre>   
	 */
	public static void testPrarmeters() {
        //1.Read in a triangle mesh from an input file with
        //  format ASCII UCD generated by Gridgen
        MeshReader reader = new MeshReader("triangle.grd");
        Mesh mesh = reader.read2DMesh();
        //Compute geometry relationship of nodes and elements
        mesh.computeNodeBelongsToElements();

        //2.Mark border types
        HashMap<NodeType, MathFunc> mapNTF =
                new HashMap<NodeType, MathFunc>();
		//Robin type on boundary x=3.0 of \Omega
		mapNTF.put(NodeType.Robin, new AbstractMathFunc("x","y"){
			@Override
			public double apply(Variable v) {
				if(3.0-v.get("x") < 0.01)
					return 1.0; //this is Robin condition
				else
					return -1.0;
			}
		});
        mapNTF.put(NodeType.Dirichlet, null);
        mesh.markBorderNode(mapNTF);

        //3.Use element library to assign degrees of
        //  freedom (DOF) to element
        ElementList eList = mesh.getElementList();
        FELinearTriangle feLT = new FELinearTriangle();
        for(int i=1;i<=eList.size();i++)
            feLT.assignTo(eList.at(i));

        //4.Weak form. We use WeakFormBuilder to define weak form
		WeakFormBuilder wfBuilder = new WeakFormBuilder() {
			/**
			 * Override this function to define weak form
			 */
			@Override
			public MathFunc makeExpression(Element e, Type type) {
				ScalarShapeFunction u = getScalarTrial();
				ScalarShapeFunction v = getScalarTest();
				MathFunc ff = getParam("f",e);
				MathFunc fk = getParam("k",e);
				MathFunc fc = getParam("c",e);
				MathFunc fd = getParam("d",e);
				MathFunc fg = getParam("g",e);
				switch(type) {
				case LHS_Domain:
					return fk.M( grad(u,"x","y").dot(grad(v,"x","y")) ).A(
						fc.M(u).M(v) );  //(k*grad(u),grad(v))_\Omega + (c*u,v)_\Omega
				case LHS_Border:
					return fd.M(u.M(v)); //(d*u,v)_\Gamma_N
				case RHS_Domain:
					return ff.M(v);      //(f,v)_\Omega
				case RHS_Border:
					return fg.M(v);      //(g,v)_\Gamma_N
				}
				return null;
			}
		};
		//-2*(x^2+y^2)+36
        wfBuilder.addParam("f",x.M(x).A(y.M(y)).M(-2.0).A(36.0));
        wfBuilder.addParam("k",C(0.1));
        wfBuilder.addParam("c",C1);
        wfBuilder.addParam("d",C1);
        wfBuilder.addParam("g",x.M(x).A(y.M(y)));
		WeakForm weakForm = wfBuilder.getScalarWeakForm();
		
        //5.Assembly process
        AssemblerScalar assembler =
                new AssemblerScalar(mesh, weakForm);
        assembler.assemble();
        Matrix stiff = assembler.getStiffnessMatrix();
        Vector load = assembler.getLoadVector();
        //Boundary condition
        assembler.imposeDirichletCondition(C0);

        //6.Solve linear system
        SolverJBLAS solver = new SolverJBLAS();
        Vector u = solver.solveDGESV(stiff, load);
        System.out.println("u=");
        for(int i=1;i<=u.getDim();i++)
            System.out.println(String.format("%.3f", u.get(i)));

        //7.Output results to an Techplot format file
        MeshWriter writer = new MeshWriter(mesh);
        writer.writeTechplot("./tutorial/UseWeakFormBuilder3.dat", u);
	}
	public static void main(String[] args) {
		triangleMesh();
		rectangleTest();
		testPrarmeters();
	}
		
}