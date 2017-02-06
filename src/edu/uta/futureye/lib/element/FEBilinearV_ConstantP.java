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
import edu.uta.futureye.core.intf.VecFiniteElement;
import edu.uta.futureye.function.FMath;
import edu.uta.futureye.function.SingleVarFunc;
import edu.uta.futureye.function.basic.FX;
import edu.uta.futureye.function.basic.SpaceVectorFunction;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.function.intf.VectorMathFunc;
import edu.uta.futureye.util.container.VertexList;
import static edu.uta.futureye.function.FMath.*;

/**
 * 2D Q1/P0 Element
 * -Continuous bilinear velocity
 * -Piecewise constant pressure
 * 
 * * Velocity: Bilinear shape function: SFBilinearLocal2D
 * * 速度：四边形局部坐标，双线性函数
 * 
 * 4----3
 * |    |
 * |    |
 * 1----2
 * 
 * NV = NV(r,s) = NV( r(x,y), s(x,y) )
 * NV1 = (1-r)*(1-s)/4
 * NV2 = (1+r)*(1-s)/4
 * NV3 = (1+r)*(1+s)/4
 * NV4 = (1-r)*(1+s)/4
 * 
 * * Pressure: Piecewise constant shape function: SFConstant1
 * * 压强：分片常数型函数
 * NP=1
 * 
 * * 2D vector valued shape functions
 * * 二维单元上的形函数，速度压强共9个自由度：
 * Ni = (u1,u2,p)', i=1,...,9
 * 
 * N1  =  (NV1, 0, 0)'
 * N2  =  (NV2, 0, 0)'
 * N3  =  (NV3, 0, 0)'
 * N4  =  (NV4, 0, 0)'
 * N5  =  (0, NV1, 0)'
 * N6  =  (0, NV2, 0)'
 * N7  =  (0, NV3, 0)'
 * N8 =   (0, NV4, 0)'
 * N9 =   (0, 0, NP)'
 *
 * @author liuyueming
 */
public class FEBilinearV_ConstantP implements VecFiniteElement {
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
	VectorMathFunc[] shapeFuncs = new VectorMathFunc[4+4+1];
	MathFunc jac;

	public FEBilinearV_ConstantP() {
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
		MathFunc[] bilinearSF = new MathFunc[4];
		bilinearSF[0] = (1-r)*(1-s)/4;
		bilinearSF[1] = (1+r)*(1-s)/4;
		bilinearSF[2] = (1+r)*(1+s)/4;
		bilinearSF[3] = (1-r)*(1+s)/4;

		//coordinate transform
		x = x1*bilinearSF[0] + x2*bilinearSF[1] + x3*bilinearSF[2] + x4*bilinearSF[3];
		y = y1*bilinearSF[0] + y2*bilinearSF[1] + y3*bilinearSF[2] + y4*bilinearSF[3];
		
		map = new HashMap<String, MathFunc>();
		map.put("x", x);
		map.put("y", y);
		// Jacobian Matrix = (r[0] r[1]) = (x_r, x_s)
		//                   (r[2] r[3])   (y_r, y_s)
		jac = x.diff("r")*y.diff("s") - y.diff("r")*x.diff("s");
		
		
		shapeFuncs[0] = new SpaceVectorFunction(bilinearSF[0], C0, C0);
		shapeFuncs[1] = new SpaceVectorFunction(bilinearSF[1], C0, C0);
		shapeFuncs[2] = new SpaceVectorFunction(bilinearSF[2], C0, C0);
		shapeFuncs[3] = new SpaceVectorFunction(bilinearSF[3], C0, C0);
		shapeFuncs[4] = new SpaceVectorFunction(C0, bilinearSF[0], C0);
		shapeFuncs[5] = new SpaceVectorFunction(C0, bilinearSF[1], C0);
		shapeFuncs[6] = new SpaceVectorFunction(C0, bilinearSF[2], C0);
		shapeFuncs[7] = new SpaceVectorFunction(C0, bilinearSF[3], C0);
		shapeFuncs[8] = new SpaceVectorFunction(C0, C0, C1);
		
		
		
	}

	@Override
	public VectorMathFunc[] getShapeFunctions() {
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
		return new FELinearLine2D();
	}

}
