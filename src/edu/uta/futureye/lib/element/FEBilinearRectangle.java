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
import edu.uta.futureye.function.SingleVarFunc;
import edu.uta.futureye.function.basic.FX;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.util.container.VertexList;

/**
How to get derivative r_x, r_y, s_x, s_y:

f(x,y) = g(r,s)
f_x = g_r*r_x + g_s*s_x  ---(1)
f_y = g_r*r_y + g_s*s_y  ---(2)

for (1), let f=x and f=y we get tow equations, solve them:
(x_r x_s)   (r_x)   (1)
(y_r y_s) * (s_x) = (0)

similarly, for (2):
(x_r x_s)   (r_y)   (0)
(y_r y_s) * (s_y) = (1)

        (x_r x_s)
Let J = (y_r y_s)

from the above four equations, we have:
 (r_x r_y)
 (s_x s_y) = inv(J)
 */
public class FEBilinearRectangle implements FiniteElement {
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
	FX x3 = new FX("x3");
	FX x4 = new FX("x4");
	FX y1 = new FX("y1");
	FX y2 = new FX("y2");
	FX y3 = new FX("y3");
	FX y4 = new FX("y4");
	RectAreaCoordR r = new RectAreaCoordR();
	RectAreaCoordS s = new RectAreaCoordS();
	//Construct a function with the coordinate of points in an element as parameters
	String[] argsOrder = new String[]{x1,x2,x3,x4,y1,y2,y3,y4,r,s};
	
	MathFunc x;
	MathFunc y;
	Map<String, MathFunc> map;
	public int nDOFs = 4;
	MathFunc[] shapeFuncs = new MathFunc[4];
	MathFunc jac;

	public FEBilinearRectangle() {
		/**
		 * shape functions
		 *  s
		 *  ^
		 *  |
		 *  |
		 * 
		 *  4-----3
		 *  |     |
		 *  |     |
		 *  1-----2  --> r
		 * -1  0  1
		 *
		 * N1 = (1-r)*(1-s)/4
		 * N2 = (1+r)*(1-s)/4
		 * N3 = (1+r)*(1+s)/4
		 * N4 = (1-r)*(1+s)/4
		 */
		shapeFuncs[0] = (1-r)*(1-s)/4;
		shapeFuncs[1] = (1+r)*(1-s)/4;
		shapeFuncs[2] = (1+r)*(1+s)/4;
		shapeFuncs[3] = (1-r)*(1+s)/4;

		//coordinate transform
		x = x1*shapeFuncs[0] + x2*shapeFuncs[1] + x3*shapeFuncs[2] + x4*shapeFuncs[3];
		y = y1*shapeFuncs[0] + y2*shapeFuncs[1] + y3*shapeFuncs[2] + y4*shapeFuncs[3];
		
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

	public void assignTo(Element e) {
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
		return new FELinearLine2D();
	}

}
