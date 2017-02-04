package edu.uta.futureye.lib.element;


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

import edu.uta.futureye.core.DOF;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Vertex;
import edu.uta.futureye.core.intf.FiniteElement;
import edu.uta.futureye.function.FMath;
import static edu.uta.futureye.function.FMath.*;
import edu.uta.futureye.function.SingleVarFunc;
import edu.uta.futureye.function.basic.FX;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.util.container.VertexList;

/**
 * Linear line element in a 1D space
 * 
 */
public class FELinearLine1D implements FiniteElement {
	public class LineCoordR extends SingleVarFunc {
		public LineCoordR() {
			super("r", "r");
		}

		@Override
		public double apply(double... args) {
			return args[this.argIdx];
		}
		
		/**
		 *  x = x1*N1 + x2*N2
		 *    = x1*(1-r)/2 + x2*(1+r)/2
		 *    = [ x1+x2 + (x2-x1)*r ]/2
		 *  =>
		 *  r = [2*x - (x1+x2)]/(x2-x1) 
		 *  r_x = 2/(x2-x1)
		 */
		@Override
		public MathFunc diff(String varName) {
			if(varName.equals("r"))
				return FMath.C1;
			if(varName.equals("x"))
				return 2.0/(x2-x1); //=1/jac
			else
				return FMath.C0;
		}

		public String toString() {
			return this.varName;
		}
		
		public String getExpr() {
			return this.varName;
		}

		@Override
		public void bytecodeGen(MethodVisitor mv, Map<String, Integer> argsMap,
				int argsStartPos, Map<MathFunc, Integer> funcRefsMap,
				String clsName) {
			mv.visitIntInsn(Opcodes.ALOAD, argsStartPos);
			Integer argIdx = argsMap.get(varName);
			if(argIdx == null) throw new RuntimeException("Index of "+varName+" is null!");
			mv.visitLdcInsn(argIdx);
			mv.visitInsn(Opcodes.DALOAD);
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

	}
	
	FX x1 = new FX("x1");
	FX x2 = new FX("x2");
	LineCoordR r = new LineCoordR();
	//Construct a function with the coordinate of points in an element as parameters
	String[] argsOrder = new String[]{x1,x2,r};
	
	MathFunc x;
	MathFunc y;
	Map<String, MathFunc> map;
	public int nDOFs = 2;
	MathFunc[] shapeFuncs = new MathFunc[nDOFs];
	MathFunc jac;

	public FELinearLine1D() {
		/**
		 * shape functions
		 * 
		 *  1-----2  -->r
		 * -1  0  1
		 * 
		 * N1 = (1-r)/2
		 * N2 = (1+r)/2
		 * 
		 * @param funID = 1,2
		 * 
		 */
		shapeFuncs[0] = 0.5*(1-r);
		shapeFuncs[1] = 0.5*(1+r);

		//coordinate transform
		x = x1*shapeFuncs[0] + x2*shapeFuncs[1];
		
		map = new HashMap<String, MathFunc>();
		map.put("x", x);
		
		/**  
		 *  Compute 1D determinant of Jacobian matrix
		 *  1D: det(Jac) = x_r
		 *  2D boundary: det(Jac)= sqrt(x_r^2 + y_r^2)
		 */
		jac = x.diff("r");
	}

	@Override
	public MathFunc[] getShapeFunctions() {
		return this.shapeFuncs;
	}

	@Override
	public int getNumberOfDOFs() {
		return this.nDOFs;
	}

	@Override
	public Map<String, MathFunc> getCoordTransMap() {
		return this.map;
	}

	@Override
	public String[] getArgsOrder() {
		return this.argsOrder;
	}
	
	@Override
	public MathFunc getJacobian() {
		return this.jac;
	}

	public void assignTo(Element e) {
		e.clearAllDOF();
		VertexList vertices = e.vertices();
		for(int j=1;j<=vertices.size();j++) {
			Vertex v = vertices.at(j);
			//Assign shape function to DOF
			DOF dof = new DOF(
						j, //Local DOF index
						v.globalNode().getIndex(), //Global DOF index, take global node index
						null //Shape function is no longer used?  
						);
			e.addNodeDOF(j, dof);
		}
	}

	@Override
	public FiniteElement getBoundaryFE() {
		return null;
	}

}
