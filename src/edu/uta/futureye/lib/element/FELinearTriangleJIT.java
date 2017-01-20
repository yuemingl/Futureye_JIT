package edu.uta.futureye.lib.element;


import java.util.HashMap;
import java.util.Map;

import org.objectweb.asm.MethodVisitor;

import com.sun.xml.internal.ws.org.objectweb.asm.Opcodes;

import edu.uta.futureye.bytecode.CompiledFunc;
import edu.uta.futureye.core.intf.FiniteElement;
import edu.uta.futureye.function.FMath;
import edu.uta.futureye.function.SingleVarFunc;
import edu.uta.futureye.function.basic.FX;
import edu.uta.futureye.function.intf.MathFunc;

public class FELinearTriangleJIT implements FiniteElement {
	public class TriAreaCoordR extends SingleVarFunc {
		public TriAreaCoordR() {
			super("r", "r");
		}

		@Override
		public double apply(double... args) {
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
		
		public String getExpr() {
			return this.varName;
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
	public class TriAreaCoordS extends SingleVarFunc {
		public TriAreaCoordS() {
			super("s", "s");
		}

		@Override
		public double apply(double... args) {
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
		
		public String getExpr() {
			return this.varName;
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

	//Construct a function with the coordinate of points in an element as parameters
	String[] argsOrder = new String[]{"x1","x2","x3","y1","y2","y3","r","s","t"};
	FX x1 = new FX("x1");
	FX x2 = new FX("x2");
	FX x3 = new FX("x3");
	FX y1 = new FX("y1");
	FX y2 = new FX("y2");
	FX y3 = new FX("y3");
	MathFunc x;
	MathFunc y;
	Map<String, MathFunc> map;
	public int nDOFs = 3;
	MathFunc[] shapeFuncs = new MathFunc[3];
	MathFunc jac;

	public FELinearTriangleJIT() {
		TriAreaCoordR r = new TriAreaCoordR();
		TriAreaCoordS s = new TriAreaCoordS();
		
		//shape functions
		shapeFuncs[0] = r;
		shapeFuncs[1] = s;
		shapeFuncs[2] = 1 - r - s;
		
		x = x1*r + x2*s + x3*(1-r-s);
		y = y1*r + y2*s + y3*(1-r-s);
		map = new HashMap<String, MathFunc>();
		map.put("x", x);
		map.put("y", y);
		// Jacobian Matrix = (r[0] r[1]) = (x_r, x_s)
		//                   (r[2] r[3])   (y_r, y_s)
		jac = x.diff("r")*y.diff("s") - y.diff("r")*x.diff("s");
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

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		FELinearTriangleJIT t = new FELinearTriangleJIT();
		TriAreaCoordR r = t.new TriAreaCoordR();
		System.out.println(r);
	}

}
