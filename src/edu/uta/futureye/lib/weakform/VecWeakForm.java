package edu.uta.futureye.lib.weakform;

import java.util.Map;

import edu.uta.futureye.bytecode.CompiledFunc;
import edu.uta.futureye.core.intf.FiniteElement;
import edu.uta.futureye.core.intf.LHSExpr;
import edu.uta.futureye.core.intf.LHSVecExpr;
import edu.uta.futureye.core.intf.RHSExpr;
import edu.uta.futureye.core.intf.RHSVecExpr;
import edu.uta.futureye.core.intf.VecFiniteElement;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.function.intf.VectorMathFunc;

public class VecWeakForm {
	VecFiniteElement fe;
//	LHSExpr lhsExpr;
//	RHSExpr rhsExpr;
	
	MathFunc jac;
	MathFunc[][] matLHS;
	MathFunc[] vecRHS;

	CompiledFunc cjac;
	CompiledFunc[][] clhs;
	CompiledFunc[] crhs;

	public VecWeakForm(VecFiniteElement fe, LHSVecExpr lhsExpr, RHSVecExpr rhsExpr) {
		this.fe = fe;
//		this.lhsExpr =  lhsExpr;
//		this.rhsExpr = rhsExpr;
		this.jac = fe.getJacobian();
 
		int nDOFs = this.fe.getNumberOfDOFs();
		VectorMathFunc[] shapeFuncs = fe.getShapeFunctions();
		Map<String, MathFunc> map = fe.getCoordTransMap();
		matLHS = new MathFunc[nDOFs][nDOFs];
		vecRHS = new MathFunc[nDOFs];

	}

	public void compile() {
		String[] argsOrder = fe.getArgsOrder();
		jac.compileToStaticField(true);
		cjac = jac.compileWithASM(argsOrder);

		int nDOFs = this.fe.getNumberOfDOFs();
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
	
	public CompiledFunc getCompiledJac() {
		return this.cjac;
	}
	
	public VecFiniteElement getFiniteElement() {
		return this.fe;
	}
}
