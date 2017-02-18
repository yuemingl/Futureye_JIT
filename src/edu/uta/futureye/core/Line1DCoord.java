/**
 * Copyright (c) 2010, nkliuyueming@gmail.com. All rights reserved.
 * 
 * 
 */
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
 * 1D line local coordinate.
 * 
 */
public class Line1DCoord implements CoordTrans {
	MathFunc x1;
	MathFunc x2;

	MathFunc r;

	MathFunc x;
	HashMap<String, MathFunc> map;

	MathFunc jac;

	public Line1DCoord(MathFunc x1, MathFunc x2) {
		this.x1 = x1;
		this.x2 = x2;

		this.r = new Line1DCoordR();
		
		MathFunc N1 = (1-r)/2;
		MathFunc N2 = (1+r)/2;

		//coordinate transform
		this.x = x1*N1 + x2*N2;

		this.map = new HashMap<String, MathFunc>();
		this.map.put("x", x);

		/**
		 *  Compute 1D determinant of Jacobian matrix
		 *  1D: det(Jac) = x_r
		 */
		jac = x.diff("r");
	}

	public MathFunc getCoordR() {
		return this.r;
	}

	@Override
	public MathFunc[] getCoords() {
		return new MathFunc[]{r};
	}

	@Override
	public MathFunc getJacobian() {
		return this.jac;
	}
	
	@Override
	public HashMap<String, MathFunc> getCoordTransMap() {
		return this.map;
	}

	public class Line1DCoordR extends SingleVarFunc {
		public Line1DCoordR() {
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
}
