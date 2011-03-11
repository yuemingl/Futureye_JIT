package edu.uta.futureye.prostate;

import edu.uta.futureye.algebra.Solver;
import edu.uta.futureye.algebra.SparseVector;
import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.function.basic.Vector2Function;
import edu.uta.futureye.lib.assembler.AssemblerScalar;
import edu.uta.futureye.lib.weakform.WeakFormDerivative;

public class Tools {
	public static Vector computeDerivative(Mesh mesh, Vector U, String varName) {
		mesh.clearBorderNodeMark();
		
		WeakFormDerivative weakForm = new WeakFormDerivative(varName);
		weakForm.setParam(new Vector2Function(U));
		
		AssemblerScalar assembler = new AssemblerScalar(mesh, weakForm);
		assembler.assemble();
		Matrix stiff = assembler.getStiffnessMatrix();
		Vector load = assembler.getLoadVector();
		
		Solver solver = new Solver();
		Vector w = solver.solve(stiff, load);
		return w;
	}
	
	public static Vector computeLaplace2D(Mesh mesh, Vector U) {
		Vector ux = Tools.computeDerivative(mesh,U,"x");
		Vector uy = Tools.computeDerivative(mesh,U,"y");
		Vector uxx = Tools.computeDerivative(mesh,ux,"x");
		Vector uyy = Tools.computeDerivative(mesh,uy,"y");
		Vector LpU = SparseVector.axpy(1.0, uxx, uyy);
		return LpU;
	}
}
