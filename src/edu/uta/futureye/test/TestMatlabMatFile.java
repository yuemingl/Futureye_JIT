package edu.uta.futureye.test;

import java.io.IOException;
import java.util.ArrayList;

import com.jmatio.io.MatFileWriter;
import com.jmatio.types.MLArray;
import com.jmatio.types.MLChar;
import com.jmatio.types.MLDouble;
import com.jmatio.types.MLSparse;

import edu.uta.futureye.algebra.SparseMatrixColMajor;
import edu.uta.futureye.algebra.SparseMatrixRowMajor;
import edu.uta.futureye.algebra.SparseVectorHashMap;
import edu.uta.futureye.algebra.intf.SparseMatrix;
import edu.uta.futureye.algebra.intf.SparseVector;
import edu.uta.futureye.io.MatlabMatFileWriter;
import edu.uta.futureye.util.Array2String;

public class TestMatlabMatFile {
	
	public static void test1() {
		 //1. First create example arrays
		 double[] src = new double[] { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 };
		 MLDouble mlDouble = new MLDouble( "double_arr", src, 3 );
		 MLChar mlChar = new MLChar( "char_arr", "I am dummy" );
		         
		 //2. write arrays to file
		 ArrayList<MLArray> list = new ArrayList<MLArray>();
		 list.add( mlDouble );
		 list.add( mlChar );
		 
		 try {
			new MatFileWriter( "mat_file.mat", list );
		} catch (IOException e) {
			e.printStackTrace();
		}		
	}
	
	public static void test2() {
		int[] dims = new int[2];
		dims[0] = 3;
		dims[1] = 3;
		MLSparse ms = new MLSparse("SM1", dims, 0, 6);
		ms.setReal(1.0, 0, 0);
		ms.setReal(2.0, 1, 0);
		ms.setReal(3.0, 1, 1);
		ms.setReal(4.0, 2, 0);
		ms.setReal(5.0, 2, 1);
		ms.setReal(6.0, 2, 2);
		System.out.println(ms.toString());
		System.out.println(Array2String.convert(ms.getIR()));
		System.out.println(Array2String.convert(ms.getJC()));

		ArrayList<MLArray> list = new ArrayList<MLArray>();
		list.add( ms );
		try {
			new MatFileWriter( "mat_file2.mat", list );
		} catch (IOException e) {
			e.printStackTrace();
		}	
	}
	
	public static void test3() {
		SparseMatrix SA = new SparseMatrixRowMajor(3,3);
		SparseMatrix SB = new SparseMatrixColMajor(3,3);
		for(int i=1;i<=3;i++) {
			for(int j=1;j<=2;j++) {
				SA.set( j,i, i+3*(j-1));
				SB.set( j,i, i+3*(j-1));
			}
		}
		SA.set(3, 1, 2.0);
		SB.set(3, 1, 2.0);
		SA.print();
		SB.print();
		
		SparseVector SV = new SparseVectorHashMap(3);
		SV.set(2, 3.0);
		
		MatlabMatFileWriter mw = new MatlabMatFileWriter();

		mw.addSparseMatrix(SA);
		mw.addSparseMatrix(SB);
		mw.addSparseVector(SV);
		
		mw.addMatrix(SA.setName("FullMatrix1"));
		mw.addMatrix(SB.setName("FullMatrix2"));
		mw.addVector(SV.setName("FullVector"));
		
		mw.writeFile("mat_file3.mat");
	}
	
	public static void main(String[] args) {
		//test1();
		//test2();
		test3();

	}
}
