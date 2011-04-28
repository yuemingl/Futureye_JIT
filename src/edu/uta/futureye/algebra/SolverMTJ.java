package edu.uta.futureye.algebra;

import java.util.Map;
import java.util.Map.Entry;

import no.uib.cipr.matrix.sparse.CG;
import no.uib.cipr.matrix.sparse.IterativeSolverNotConvergedException;
import edu.uta.futureye.algebra.intf.Vector;

public class SolverMTJ {
	
	/**
	 * Matrix-Toolkits-Java(MTJ) interface
	 * 
	 */
	public Vector solveCG(SparseMatrix A, SparseVector b, 
			SparseVector x) {
		int dim = b.getDim();
		no.uib.cipr.matrix.sparse.SparseVector template = 
			new no.uib.cipr.matrix.sparse.SparseVector(dim);
		CG cg = new CG(template);
		
		no.uib.cipr.matrix.sparse.CompRowMatrix A2 = 
			new no.uib.cipr.matrix.sparse.CompRowMatrix(
					A.getRowDim(),A.getColDim(),A.getColIndex());
		
		Map<Integer, Map<Integer, Double>> data = A.getAll();
		for(Entry<Integer,Map<Integer,Double>> row : data.entrySet()) {
			int nRow = row.getKey();
			for(Entry<Integer,Double> col : row.getValue().entrySet()) {
				int nCol = col.getKey();
				A2.set(nRow-1, nCol-1, col.getValue());
			}
			row.getValue().clear();
		}
		data.clear();
		
		no.uib.cipr.matrix.sparse.SparseVector b2 = 
			new no.uib.cipr.matrix.sparse.SparseVector(dim);
		no.uib.cipr.matrix.sparse.SparseVector x2 = 
			new no.uib.cipr.matrix.sparse.SparseVector(dim);
		Map<Integer,Double> bData = b.getAll();
		for(Entry<Integer,Double> ety : bData.entrySet()) {
			b2.set(ety.getKey()-1,ety.getValue());
		}
		for(int i=0;i<dim;i++) {
			x2.set(i, 0.01);
		}
		
		try {
			long begin = System.currentTimeMillis();
			cg.solve(A2, b2, x2);
			long end = System.currentTimeMillis();
			System.out.println("Solve time used:"+(end-begin));
			System.out.println("Iterations: "+cg.getIterationMonitor().iterations());

		} catch (IterativeSolverNotConvergedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		
		for(int i=1;i<=dim;i++) {
			x.set(i, x2.get(i-1));
		}
		
		return x;
    }
}
