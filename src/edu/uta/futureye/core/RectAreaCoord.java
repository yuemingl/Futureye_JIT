package edu.uta.futureye.core;

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

import edu.uta.futureye.core.intf.CoordTrans;
import edu.uta.futureye.function.FMath;
import edu.uta.futureye.function.SingleVarFunc;
import edu.uta.futureye.function.intf.MathFunc;

/**
 * Rectangle area coordinate r,s (bilinear).
 *
 * Remark: 
 * How to get the derivative: r_x, r_y, s_x, s_y:
 * 
 * f(x,y) = g(r,s)
 * f_x = g_r*r_x + g_s*s_x  ---(1)
 * f_y = g_r*r_y + g_s*s_y  ---(2)
 * 
 * for (1), let f=x and f=y we get tow equations, solve them:
 * (x_r x_s)   (r_x)   (1)
 * (y_r y_s) * (s_x) = (0)
 *
 * similarly, for (2):
 * (x_r x_s)   (r_y)   (0)
 * (y_r y_s) * (s_y) = (1)
 * 
 * Let J = (x_r x_s)
 *         (y_r y_s)
 * 
 * from the above four equations, we have:
 *  (r_x r_y) = inv(J)
 *  (s_x s_y)
 *
 */
public class RectAreaCoord implements CoordTrans {
	MathFunc x1;
	MathFunc x2;
	MathFunc x3;
	MathFunc x4;
	MathFunc y1;
	MathFunc y2;
	MathFunc y3;
	MathFunc y4;

	RectAreaCoordR r;
	RectAreaCoordS s;

	MathFunc x;
	MathFunc y;
	HashMap<String, MathFunc> map;

	MathFunc jac;

	/**
	 * 
	 * @param x1
	 * @param x2
	 * @param x3
	 * @param x4
	 * @param y1
	 * @param y2
	 * @param y3
	 * @param y4
	 */
	public RectAreaCoord(MathFunc x1, MathFunc x2, MathFunc x3, MathFunc x4,
			MathFunc y1, MathFunc y2, MathFunc y3, MathFunc y4) {
		this.x1 = x1;
		this.x2 = x2;
		this.x3 = x3;
		this.x4 = x4;
		this.y1 = y1;
		this.y2 = y2;
		this.y3 = y3;
		this.y4 = y4;

		this.r = new RectAreaCoordR();
		this.s = new RectAreaCoordS();

		MathFunc N1 = (1-r)*(1-s)/4;
		MathFunc N2 = (1+r)*(1-s)/4;
		MathFunc N3 = (1+r)*(1+s)/4;
		MathFunc N4 = (1-r)*(1+s)/4;

		//coordinate transform
		this.x = x1*N1 + x2*N2 + x3*N3 + x4*N4;
		this.y = y1*N1 + y2*N2 + y3*N3 + y4*N4;

		this.map = new HashMap<String, MathFunc>();
		this.map.put("x", x);
		this.map.put("y", y);
		
		// Jacobian Matrix = (r[0] r[1]) = (x_r, x_s)
		//                   (r[2] r[3])   (y_r, y_s)
		this.jac = x.diff("r")*y.diff("s") - y.diff("r")*x.diff("s");
	}

	public MathFunc getCoordR() {
		return this.r;
	}

	public MathFunc getCoordS() {
		return this.s;
	}

	@Override
	public MathFunc[] getCoords() {
		return new MathFunc[]{r, s};
	}

	@Override
	public MathFunc getJacobian() {
		return this.jac;
	}
	
	@Override
	public HashMap<String, MathFunc> getCoordTransMap() {
		return this.map;
	}

	public class RectAreaCoordR extends SingleVarFunc {
		public RectAreaCoordR() {
			super("r", "r");
		}

		@Override
		public double apply(double... args) {
			return args[this.argIdx];
		}
		
		@Override
		public MathFunc diff(String varName) {
			if(varName.equals("r"))
				return FMath.C1;
			if(varName.equals("x"))
				return y.diff("s")/jac;
			else if(varName.equals("y"))
				return -x.diff("s")/jac;
			else
				return FMath.C0;
		}

		public String getExpr() {
			return this.varName;
		}

		public String toString() {
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
	public class RectAreaCoordS extends SingleVarFunc {
		public RectAreaCoordS() {
			super("s", "s");
		}

		@Override
		public double apply(double... args) {
			return args[this.argIdx];
		}

		@Override
		public MathFunc diff(String varName) {
			if(varName.equals("s"))
				return FMath.C1;
			if(varName.equals("x"))
				return -y.diff("r")/jac;
			else if(varName.equals("y"))
				return x.diff("r")/jac;
			else
				return FMath.C0;
		}

		public String getExpr() {
			return this.varName;
		}

		public String toString() {
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
}
