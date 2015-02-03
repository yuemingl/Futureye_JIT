package edu.uta.futureye.io;

import java.io.FileNotFoundException;
import java.io.IOException;

import com.jmatio.io.MatFileReader;
import com.jmatio.types.MLArray;
import com.jmatio.types.MLDouble;
import com.jmatio.types.MLInt8;
import com.jmatio.types.MLUInt8;

import edu.uta.futureye.algebra.SpaceVector;
import edu.uta.futureye.algebra.SparseMatrixRowMajor;
import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.util.FutureyeException;

public class MatlabMatFileReader {
	MatFileReader rd = null;
	
	public MatlabMatFileReader(String fileName) {
		try {
			rd = new MatFileReader(fileName);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public Vector getVector(String vecName) {
		Vector rlt = null;
		MLArray a = rd.getMLArray(vecName);
		if(a instanceof MLDouble) {
			MLDouble ad = (MLDouble)a;
			if(ad.getM()==1 || ad.getN()==1) {
				double[][] d = ad.getArray();
				if(d.length == 1) {//行向量
					rlt = new SpaceVector(d[0].length);
					for(int i=0;i<d[0].length;i++) {
						rlt.set(i+1, d[0][i]);
					}
				} else {//列向量
					rlt = new SpaceVector(d.length);
					for(int i=0;i<d.length;i++) {
						rlt.set(i+1, d[i][0]);
					}
				}
			} else {
				throw new FutureyeException(vecName+" is not a vector!");
			}
		} else {
			throw new FutureyeException("Unsported array type, name:"+vecName+" array:"+a);
		}
		return rlt;
	}
	
	public Matrix getMatrix(String matName) {
		Matrix rlt = null;
		MLArray a = rd.getMLArray(matName);
		if(a instanceof MLDouble) {
			MLDouble ad = (MLDouble)a;
			int m = ad.getM();
			int n = ad.getN();
			rlt = new SparseMatrixRowMajor(m,n);
			double[][] d = ad.getArray();
			for(int i=0;i<m;i++) {
				for(int j=0;j<n;j++) {
					rlt.set(i+1, j+1, d[i][j]);
				}
			}
		} else if(a instanceof MLUInt8) {
			MLUInt8 ad = (MLUInt8)a;
			int m = ad.getM();
			int n = ad.getN();
			rlt = new SparseMatrixRowMajor(m,n);
			byte[][] d = ad.getArray();
			for(int i=0;i<m;i++) {
				for(int j=0;j<n;j++) {
					rlt.set(i+1, j+1, d[i][j]);
				}
			}		
		} else {
			throw new FutureyeException("Unsported array type, name:"+matName+" array:"+a);
		}
		return rlt;
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		MatlabMatFileReader r = new MatlabMatFileReader("./HumanReal/interplate_data.mat");
		Vector v = r.getVector("XI");
		Matrix m = r.getMatrix("VVI");
		v.print();
		m.print();
	}

}
