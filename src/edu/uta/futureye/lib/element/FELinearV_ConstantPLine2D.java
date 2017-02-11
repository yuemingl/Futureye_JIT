package edu.uta.futureye.lib.element;

import static edu.uta.futureye.function.FMath.C0;
import static edu.uta.futureye.function.FMath.C1;
import static edu.uta.futureye.function.FMath.pow;
import static edu.uta.futureye.function.FMath.sqrt;

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

import edu.uta.futureye.core.intf.VecFiniteElement;
import edu.uta.futureye.function.FMath;
import edu.uta.futureye.function.SingleVarFunc;
import edu.uta.futureye.function.basic.FX;
import edu.uta.futureye.function.basic.SpaceVectorFunction;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.function.intf.VectorMathFunc;

public class FELinearV_ConstantPLine2D implements VecFiniteElement {
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
			if(varName.equals("y"))
				return 2.0/(y2-y1); //=1/jac
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
	FX y1 = new FX("y1");
	FX y2 = new FX("y2");
	LineCoordR r = new LineCoordR();
	//Construct a function with the coordinate of points in an element as parameters
	String[] argsOrder = new String[]{x1,x2,y1,y2,r};
	
	MathFunc x;
	MathFunc y;
	Map<String, MathFunc> map;
	public int nDOFs = 2+2+1;
	VectorMathFunc[] shapeFuncs = new VectorMathFunc[nDOFs];
	MathFunc jac;

	public FELinearV_ConstantPLine2D() {
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
		//shape function for velocity component
		MathFunc NV1 = 0.5*(1-r);
		MathFunc NV2 = 0.5*(1+r);

		//coordinate transform
		x = x1*NV1 + x2*NV2;
		y = y1*NV1 + y2*NV2;
		
		map = new HashMap<String, MathFunc>();
		map.put("x", x);
		
		/**  
		 *  Compute 1D determinant of Jacobian matrix
		 *  1D: det(Jac) = x_r
		 *  2D boundary: det(Jac)= sqrt(x_r^2 + y_r^2)
		 */
		jac = sqrt(pow(x.diff("r"),2) + pow(y.diff("r"),2));

		shapeFuncs[0] = new SpaceVectorFunction(NV1, C0, C0);
		shapeFuncs[1] = new SpaceVectorFunction(NV2, C0, C0);
		shapeFuncs[2] = new SpaceVectorFunction(C0, NV1, C0);
		shapeFuncs[3] = new SpaceVectorFunction(C0, NV2, C0);
		shapeFuncs[4] = new SpaceVectorFunction(C0, C0, C1);

	}

	@Override
	public int getNumberOfDOFs() {
		return this.nDOFs;
	}

	@Override
	public VectorMathFunc[] getShapeFunctions() {
		return this.shapeFuncs;
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

	@Override
	public VecFiniteElement getBoundaryFE() {
		return null;
	}

	@Override
	public boolean isDOFCoupled(int idx1, int idx2) {
		if(idx1 == 4 || idx2 == 4)
			return true;
		if(idx1 <= 1 && idx2 >=2)
			return false;
		if(idx2 <= 1 && idx1 >=2)
			return false;
		return true;
	}

}
