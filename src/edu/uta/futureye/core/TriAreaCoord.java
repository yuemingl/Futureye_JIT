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

import edu.uta.futureye.function.FMath;
import edu.uta.futureye.function.SingleVarFunc;
import edu.uta.futureye.function.intf.MathFunc;

/**
 * Triangular area coordinate
 * r,s,t
 * t=1-r-s
 * @author yueming.liu
 *
 */
public class TriAreaCoord {
	MathFunc x1;
	MathFunc x2;
	MathFunc x3;
	MathFunc y1;
	MathFunc y2;
	MathFunc y3;
	
	TriAreaCoordR r;
	TriAreaCoordS s;
	
	MathFunc x;
	MathFunc y;
	HashMap<String, MathFunc> map;
	
	MathFunc jac;

	/**
	 * 
	 * @param x1
	 * @param x2
	 * @param x3
	 * @param y1
	 * @param y2
	 * @param y3
	 */
	public TriAreaCoord(MathFunc x1, MathFunc x2, MathFunc x3, MathFunc y1, MathFunc y2, MathFunc y3) {
		this.x1 = x1;
		this.x2 = x2;
		this.x3 = x3;
		this.y1 = y1;
		this.y2 = y2;
		this.y3 = y3;
		
		this.r = new TriAreaCoordR();
		this.s = new TriAreaCoordS();
		
		this.x = x1*r + x2*s + x3*(1-r-s);
		this.y = y1*r + y2*s + y3*(1-r-s);
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
	
	/**
	 * 
	 * @return 1-r-s
	 */
	public MathFunc getCoordT() {
		return 1 - r - s;
	}
	
	public MathFunc getJacobian() {
		return this.jac;
	}
	
	public HashMap<String, MathFunc> getCoordTransMap() {
		return this.map;
	}
	
	public class TriAreaCoordR extends SingleVarFunc {
		public TriAreaCoordR() {
			super("r", "r");
		}
	
		@Override
		public double apply(double... args) {
			//this.argIdx is wrong if we don't define BCEL bytecodeGen
			//
			return args[this.argIdx];
		}
		
		@Override
		public MathFunc diff(String varName) {
			if(varName.equals("r"))
				return FMath.C1;
			if(varName.equals("x")) //r_x
				return (y2-y3)/jac;
			else if(varName.equals("y")) //r_y
				return (x3-x2)/jac;
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
			mv.visitLdcInsn(argsMap.get(varName));
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
	public class TriAreaCoordS extends SingleVarFunc {
		public TriAreaCoordS() {
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
			if(varName.equals("x")) //s_x
				return (y3-y1)/jac;
			else if(varName.equals("y")) //s_y
				return (x1-x3)/jac;
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
			mv.visitLdcInsn(argsMap.get(varName));
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
