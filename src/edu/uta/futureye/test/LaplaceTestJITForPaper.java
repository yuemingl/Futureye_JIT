package edu.uta.futureye.test;

import static edu.uta.futureye.function.FMath.C0;
import static edu.uta.futureye.function.FMath.grad;
import static edu.uta.futureye.function.FMath.r;
import static edu.uta.futureye.function.FMath.s;
import static edu.uta.futureye.function.FMath.x;
import static edu.uta.futureye.function.FMath.y;

import java.util.HashMap;
import java.util.Map;

import org.apache.bcel.generic.ALOAD;
import org.apache.bcel.generic.ConstantPoolGen;
import org.apache.bcel.generic.DALOAD;
import org.apache.bcel.generic.InstructionFactory;
import org.apache.bcel.generic.InstructionHandle;
import org.apache.bcel.generic.InstructionList;
import org.apache.bcel.generic.MethodGen;
import org.apache.bcel.generic.PUSH;
import org.objectweb.asm.MethodVisitor;

import com.sun.xml.internal.ws.org.objectweb.asm.Opcodes;

import edu.uta.futureye.algebra.SparseMatrixRowMajor;
import edu.uta.futureye.algebra.SparseVectorHashMap;
import edu.uta.futureye.algebra.intf.SparseMatrix;
import edu.uta.futureye.algebra.intf.SparseVector;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.algebra.solver.external.SolverJBLAS;
import edu.uta.futureye.bytecode.CompiledFunc;
import edu.uta.futureye.core.DOF;
import edu.uta.futureye.core.DOFOrder;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.NodeType;
import edu.uta.futureye.function.FMath;
import edu.uta.futureye.function.SingleVarFunc;
import edu.uta.futureye.function.basic.FX;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.function.operator.FOIntegrate;
import edu.uta.futureye.io.MeshReader;
import edu.uta.futureye.io.MeshWriter;
import edu.uta.futureye.lib.element.FELinearTriangleOld;
import edu.uta.futureye.util.MeshGenerator;
import edu.uta.futureye.util.Utils;
import edu.uta.futureye.util.container.DOFList;
import edu.uta.futureye.util.container.ElementList;


/**
 * compileWithASM and benchmark for paper
 *
 * <blockquote><pre>
 * Problem:
 *   -\Delta{u} = f
 *   u(x,y)=0, (x,y) \in \partial{\Omega}
 * where
 *   \Omega = [-3,3]*[-3,3]
 *   f = 2 * pi * pi * sin ( pi * x ) * sin ( pi * y )
 * Solution:
 *   u = (x^2-9)*(y^2-9)
 * </blockquote></pre>
 * 
 * @author liuyueming
 */

public class LaplaceTestJITForPaper {
	
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
			//argIdx is not correct??????
			//if we do not override 'InstructionHandle bytecodeGen(..)'
			//this function will be call and argIdx need to be figured out
			//this.argIdx is set in 'BytecodeUtils.genClass()'
			return args[this.argIdx];
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
		@Override
		public MathFunc setArgIdx(Map<String, Integer> argsMap) {
			this.argIdx = argsMap.get(varName);
			return this;
		}

		@Override
		public InstructionHandle bytecodeGen(String clsName, MethodGen mg, 
				ConstantPoolGen cp, InstructionFactory factory, 
				InstructionList il, Map<String, Integer> argsMap, 
				int argsStartPos, Map<MathFunc, Integer> funcRefsMap) {
			il.append(new ALOAD(argsStartPos));
			il.append(new PUSH(cp, argsMap.get(this.getName())));
			return il.append(new DALOAD());
		}

		@Override
		public void bytecodeGen(MethodVisitor mv, Map<String, Integer> argsMap,
				int argsStartPos, Map<MathFunc, Integer> funcRefsMap,
				String clsName) {
			mv.visitIntInsn(Opcodes.ALOAD, argsStartPos);
			mv.visitLdcInsn(argsMap.get(varName));
			mv.visitInsn(Opcodes.DALOAD);
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
			//argIdx is not correct??????
			//if we do not override 'InstructionHandle bytecodeGen(..)'
			//this function will be call and argIdx need to be figured out
			//this.argIdx is set in 'BytecodeUtils.genClass()'
			return args[this.argIdx];
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
		
		@Override
		public InstructionHandle bytecodeGen(String clsName, MethodGen mg, 
				ConstantPoolGen cp, InstructionFactory factory, 
				InstructionList il, Map<String, Integer> argsMap, 
				int argsStartPos, Map<MathFunc, Integer> funcRefsMap) {
			il.append(new ALOAD(argsStartPos));
			il.append(new PUSH(cp, argsMap.get(this.getName())));
			return il.append(new DALOAD());
		}
		
		@Override
		public void bytecodeGen(MethodVisitor mv, Map<String, Integer> argsMap,
				int argsStartPos, Map<MathFunc, Integer> funcRefsMap,
				String clsName) {
			mv.visitIntInsn(Opcodes.ALOAD, argsStartPos);
			mv.visitLdcInsn(argsMap.get(varName));
			mv.visitInsn(Opcodes.DALOAD);
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
			// 2D JacMat = (r[0] r[1]) = (x_r, x_s)
			//             (r[2] r[3])   (y_r, y_s)
			//jac changes with element, define the expression for jac with linear element
			jac = fx.diff("r")*fy.diff("s") - fy.diff("r")*fx.diff("s");
			jac.compileToStaticField(true);
			cjac = jac.compileWithASM(argsOrder);
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
//				System.out.println(">>>"+sf[j]);
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
		CompiledFunc cjac;
		
		public void compileWeakForm() {
			clhs = new CompiledFunc[nDOFs][nDOFs];
			crhs = new CompiledFunc[nDOFs];
			for(int j=0; j<nDOFs; j++) {
				for(int i=0; i<nDOFs; i++) {
					clhs[j][i] = matLHS[j][i].compileWithASM(argsOrder);
					//clhs[j][i] = matLHS[j][i].compile(argsOrder);
				}
				crhs[j] = vecRHS[j].compileWithASM(argsOrder);
				//crhs[j] = vecRHS[j].compile(argsOrder);
			}
		}
		
		public CompiledFunc[][] getCompiledLHS() {
			return clhs;
		}
		
		public CompiledFunc[] getCompiledRHS() {
			return crhs;
		}
		
		public CompiledFunc getJac() {
			return this.cjac;
		}
		
	}
	
	public void run(int nNodes) {
		int n = 51;
		boolean solveSystem = false;

		// 1.Generate mesh
		Mesh mesh = null;
		if (solveSystem) {
			MeshReader reader = new MeshReader("triangle.grd");
			mesh = reader.read2DMesh();
		} else {
			mesh = MeshGenerator.rectangle(-3, 3, -3, 3, n, n);
		}

		// Compute geometry relationship between nodes and elements
		mesh.computeNodeBelongsToElements();

		// 2.Mark border types
		HashMap<NodeType, MathFunc> mapNTF = new HashMap<NodeType, MathFunc>();
		mapNTF.put(NodeType.Dirichlet, null);
		mesh.markBorderNode(mapNTF);

		// 3.Use element library to assign degrees of
		// freedom (DOF) to element
		ElementList eList = mesh.getElementList();
		FELinearTriangleOld feLT = new FELinearTriangleOld();
		for (int i = 1; i <= eList.size(); i++)
			feLT.assignTo(eList.at(i));

		// Construct a function with the coordinate of points in an element as
		// parameters
		String[] argsOrder = new String[] { "x1", "x2", "x3", "y1", "y2", "y3",
				"r", "s", "t" };
		FELinearTriangleT fet = new FELinearTriangleT();

		// Right hand side(RHS):
		// MathFunc f = 2.0 * PI * PI * sin ( PI * x ) * sin ( PI * y );
		final MathFunc f = -2 * (x * x + y * y) + 36;

		// fet.makeWeakForm(
		// (u,v) -> grad(u,"x","y").dot(grad(v,"x","y")),
		// v -> f*v
		// );
		fet.makeWeakForm(new LHSExpr() {
			public MathFunc apply(MathFunc u, MathFunc v) {
				return grad(u, "x", "y").dot(grad(v, "x", "y"));
			}
		}, new RHSExpr() {
			public MathFunc apply(MathFunc v) {
				return f * v;
			}
		});

		long startCompile = System.currentTimeMillis();
		fet.compileWeakForm();
		System.out.println("Compile time: "
				+ (System.currentTimeMillis() - startCompile));

		CompiledFunc[][] clhs = fet.getCompiledLHS();
		CompiledFunc[] crhs = fet.getCompiledRHS();
		int nDOFs = fet.nDOFs;

		// 5.Assembly process
		double[][] A = new double[nDOFs][nDOFs];
		double[] b = new double[nDOFs];
		double[] params = new double[argsOrder.length];
		int dim = mesh.getNodeList().size();
		SparseMatrix stiff = new SparseMatrixRowMajor(dim, dim);
		SparseVector load = new SparseVectorHashMap(dim);

		long start = System.currentTimeMillis();
		int NN = nNodes / ((n - 1) * (n - 1)); // 10000*512/eList.size();
		if (solveSystem)
			NN = 1;
		for (int ii = 0; ii < NN; ii++) {
			for (Element e : eList) {
				// e.adjustVerticeToCounterClockwise();

				DOFList DOFs = e.getAllDOFList(DOFOrder.NEFV);
				double[] coords = e.getNodeCoords();
				System.arraycopy(coords, 0, params, 0, coords.length);

				// put it in intOnTriangleRefElement?
				fet.getJac().apply(params);

				for (int j = 0; j < nDOFs; j++) {
					for (int i = 0; i < nDOFs; i++) {
						A[j][i] = FOIntegrate.intOnTriangleRefElement(
								clhs[j][i], params, coords.length, 2);// 2=80.839
																		// 3=80.966,
																		// 4=80.967
					}
					b[j] = FOIntegrate.intOnTriangleRefElement(crhs[j], params,
							coords.length, 2);
				}

				if (solveSystem) {
					for (int j = 0; j < nDOFs; j++) {
						DOF dofI = DOFs.at(j + 1);
						int nGlobalRow = dofI.getGlobalIndex();
						for (int i = 0; i < nDOFs; i++) {
							DOF dofJ = DOFs.at(i + 1);
							int nGlobalCol = dofJ.getGlobalIndex();
							stiff.add(nGlobalRow, nGlobalCol, A[j][i]);
						}
						// Local load vector
						load.add(nGlobalRow, b[j]);
					}
				}
			}
		}
		System.out.println("Nodes=" + nNodes + ", Aassembly time: "
				+ (System.currentTimeMillis() - start) + "ms");

		if (solveSystem) {
			// Boundary condition
			Utils.imposeDirichletCondition(stiff, load, mesh, C0);

			// 6.Solve linear system
			SolverJBLAS solver = new SolverJBLAS();
			Vector u = solver.solveDGESV(stiff, load);
			System.out.println("u=");
			for (int i = 1; i <= u.getDim(); i++)
				System.out.println(String.format("%.3f ", u.get(i)));

			// 7.Output results to an Techplot format file
			MeshWriter writer = new MeshWriter(mesh);
			writer.writeTechplot("./tutorial/Laplace2D.dat", u);

			this.mesh = mesh;
			this.u = u;
		}
	}

    public static void main(String[] args) {
    	LaplaceTestJITForPaper ex1 = new LaplaceTestJITForPaper();
    	ex1.run(10000);
    	ex1.run(100000);
    	ex1.run(10000*512);
    	ex1.run(1000000);
    	ex1.run(10000000);
    	ex1.run(100000000);
    	/**
    	 * Nodes=10000, Aassembly time: 120ms
Nodes=100000, Aassembly time: 249ms
Nodes=1000000, Aassembly time: 1719ms
Nodes=10000000, Aassembly time: 15709ms
Nodes=100000000, Aassembly time: 147688ms


Compile time: 130
Nodes=10000, Aassembly time: 118ms
Compile time: 23
Nodes=100000, Aassembly time: 252ms
Compile time: 27
Nodes=1000000, Aassembly time: 1625ms
Compile time: 25
Nodes=10000000, Aassembly time: 15896ms
Compile time: 18
Nodes=100000000, Aassembly time: 142159ms


    	 */
    }
}
