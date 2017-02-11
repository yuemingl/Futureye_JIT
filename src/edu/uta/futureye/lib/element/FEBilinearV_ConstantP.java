package edu.uta.futureye.lib.element;


import static edu.uta.futureye.function.FMath.C0;
import static edu.uta.futureye.function.FMath.C1;

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
	public int nDOFs = 4+4+1;
	VectorMathFunc[] shapeFuncs = new VectorMathFunc[nDOFs];
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
		MathFunc N1 = (1-r)*(1-s)/4;
		MathFunc N2 = (1+r)*(1-s)/4;
		MathFunc N3 = (1+r)*(1+s)/4;
		MathFunc N4 = (1-r)*(1+s)/4;

		//coordinate transform
		x = x1*N1 + x2*N2 + x3*N3 + x4*N4;
		y = y1*N1 + y2*N2 + y3*N3 + y4*N4;
		
		map = new HashMap<String, MathFunc>();
		map.put("x", x);
		map.put("y", y);
		// Jacobian Matrix = (r[0] r[1]) = (x_r, x_s)
		//                   (r[2] r[3])   (y_r, y_s)
		jac = x.diff("r")*y.diff("s") - y.diff("r")*x.diff("s");
		
		shapeFuncs[0] = new SpaceVectorFunction(N1, C0, C0);
		shapeFuncs[1] = new SpaceVectorFunction(N2, C0, C0);
		shapeFuncs[2] = new SpaceVectorFunction(N3, C0, C0);
		shapeFuncs[3] = new SpaceVectorFunction(N4, C0, C0);
		shapeFuncs[4] = new SpaceVectorFunction(C0, N1, C0);
		shapeFuncs[5] = new SpaceVectorFunction(C0, N2, C0);
		shapeFuncs[6] = new SpaceVectorFunction(C0, N3, C0);
		shapeFuncs[7] = new SpaceVectorFunction(C0, N4, C0);
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

//	//???
//	// we need total number of nodes in the mesh to assign global index for u[2]
	
	//DOF contains local-global index 
	// no shape functions now.
	
//	public void assignTo(Element e) {
//		e.clearAllDOF();
//		if(nTotalNodes == -1 || nDOF_p == -1) {
//			FutureyeException ex = new FutureyeException("Call initDOFIndex() first!");
//			ex.printStackTrace();
//			System.exit(-1);
//		}
//		//单元结点数
//		int nNode = e.nodes.size();
//		//Assign shape function to DOF
//		for(int j=1;j<=nNode;j++) {
//			//Asign shape function to DOF
//			DOF dof_u1 = new DOF(
//					j,//Local DOF index
//					//Global DOF index, take global node index
//					e.nodes.at(j).globalIndex,
//					null
//					         );
//			dof_u1.setVVFComponent(1);
//			DOF dof_u2 = new DOF(
//					nNode+j,//Local DOF index
//					//Global DOF index, take this.nTotalNodes + global node index
//					this.nTotalNodes+e.nodes.at(j).globalIndex,
//					null
//					         );
//			dof_u2.setVVFComponent(2);
//			e.addNodeDOF(j, dof_u1);
//			//e.addNodeDOF(j, dof_u2); //???bug???
//			e.addNodeDOF(nNode+j, dof_u2);
//		}
//		
//		//Assign shape function to DOF
//		DOF dof = new DOF(
//					2*nNode+1, //Local DOF index
//					//this.nTotalNodes*2+nDOF_p, //Global DOF index for Pressure
//					this.nTotalNodes*2+this.nDOF_p, //Global DOF index for Pressure
//					shapeFun[2*nNode] //Shape function 
//					);
//		this.nDOF_p++;
//		dof.setVVFComponent(3);	
//		e.addVolumeDOF(dof);
//	}

	@Override
	public VecFiniteElement getBoundaryFE() {
		return new FELinearV_ConstantPLine2D();
	}

	@Override
	public boolean isDOFCoupled(int idx1, int idx2) {
		if(idx1 == 8 || idx2 == 8)
			return true;
		if(idx1 <= 3 && idx2 >= 4)
			return false;
		if(idx2 <= 3 && idx1 >= 4)
			return false;
		return true;
	}
}
