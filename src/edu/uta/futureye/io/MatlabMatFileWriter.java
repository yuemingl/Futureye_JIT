package edu.uta.futureye.io;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Map;
import java.util.Map.Entry;

import com.jmatio.io.MatFileWriter;
import com.jmatio.types.MLArray;
import com.jmatio.types.MLDouble;
import com.jmatio.types.MLSparse;

import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.MatrixEntry;
import edu.uta.futureye.algebra.intf.SparseMatrix;
import edu.uta.futureye.algebra.intf.SparseVector;
import edu.uta.futureye.algebra.intf.Vector;

/**
 * An interface to Matlab mat file format using JMatIO package.
 * <p>
 * With this class, matrices and vectors in FuturEye can be stored as mat file format of Matlab.
 * Use command <tt>load(fileName)</tt> in matlab to load these objects into matlab workspace for further processes.
 * 
 * @author liuyueming
 *
 */
public class MatlabMatFileWriter {
	protected ArrayList<MLArray> list = new ArrayList<MLArray>();
	
	public MatlabMatFileWriter() {
	}
	
	/**
	 * Add a matrix to the writer, the matrix will be a full matrix in matlab workspace
	 * 
	 * @param mat
	 */
	public void addMatrix(Matrix mat) {
		int[] dims = new int[2];
		dims[0] = mat.getRowDim();
		dims[1] = mat.getColDim();
		MLDouble  md = new MLDouble(mat.getName(),dims);
		for(int i=1;i<=dims[0];i++)
			for(int j=1;j<=dims[1];j++) 
				md.setReal(mat.get(i, j), i-1, j-1);
		list.add(md);
	}
	
	/**
	 * 
	 * Add a vector to the writer, the vector will be a full row-vector in matlab workspace
	 * 
	 * @param vec
	 */
	public void addVector(Vector vec) {
		int dim = vec.getDim();
		double[] data = new double[dim];
		for(int i=dim+1; --i>0;)
			data[i-1] = vec.get(i);
		MLDouble  md = new MLDouble(vec.getName(),data,1);
		list.add(md);
	}
	
	/**
	 * Add a sparse matrix to the writer, the matrix will also be a sparse matrix in matlab workspace
	 * 
	 * @param mat
	 */
	public void addSparseMatrix(SparseMatrix mat) {
		int[] dims = new int[2];
		dims[0] = mat.getRowDim();
		dims[1] = mat.getColDim();
		MLSparse sparse = new MLSparse(mat.getName(),dims,0,mat.getNonZeroNumber());
		for(MatrixEntry e : mat) {
			sparse.setReal(e.getValue(), e.getRow()-1, e.getCol()-1);
		}
		list.add(sparse);
	}
	
	/**
	 * Add a sparse vector to the writer, the vector will also be a sparse row-vector in matlab workspace
	 * 
	 * @param mat
	 */
	public void addSparseVector(SparseVector vec) {
		Map<Integer, Double> data = vec.getAll();
		int[] dims = new int[2];
		dims[0] = 1;
		dims[1] = vec.getDim();
		MLSparse sparse = new MLSparse(vec.getName(),dims,0,vec.getNonZeroNumber());
		for(Entry<Integer, Double> e : data.entrySet()) {
			int nCol = e.getKey();
			double val = e.getValue();
			sparse.setReal(val, nCol-1, nCol-1);
		}
		list.add(sparse);
	}
	
	/**
	 * Write all the matrices and vectors to the *.mat file named by <tt>fineName</tt>
	 * 
	 * @param fileName Matlab *.mat File name. If <tt>fileName</tt> is not ended with '.mat' a postfix '.mat' will be added.
	 */
	public void writeFile(String fileName) {
		try {
			if(!fileName.endsWith(".mat"))
				fileName = fileName + ".mat";
			new MatFileWriter(fileName, list);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
}
