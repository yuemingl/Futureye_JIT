package edu.uta.futureye.application;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;

import edu.uta.futureye.algebra.SpaceVector;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.util.FutureyeException;

public class DataReader {
	public static boolean debug = false;
	
	public static Vector readVector(String file) {
		return readVector(file, 3);
	}
	
	public static Vector readVector(String file, int nCol) {
		FileInputStream in;
		Vector v = null;
		try {
			in = new FileInputStream(file);

			InputStreamReader reader = new InputStreamReader(in);
			BufferedReader br = new BufferedReader(reader);
	
			String str = null;
			int N = 0;
			int dim = 1;
			while((str = br.readLine()) != null){
				if(debug)
					System.out.println(str);
				if(N == 0) {
					int n = str.indexOf("N=");
					int e = str.indexOf("E=");
					if(n != -1) {
						N = Integer.parseInt(str.substring(n+2,e-1));
						v = new SpaceVector(N);
					}
				} else {
					String[] line = str.split("(\\s)+");
					v.set(dim,Double.valueOf(line[nCol-1]));
					dim++;
					N--;
					if(N==0) break;
				}
			}
			if(N != 0) throw new FutureyeException("The dimension of the vector is wrong! file:"+file);
			//v.setDim(dim-1);
			br.close();
			in.close();

			return v;
		
		} catch (Exception e) {
			e.printStackTrace();
		}
		return null;
		
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		Vector v = readVector(".\\MouseHead\\Results\\760nm\\B12L6\\alpha_omega_BL01_ext_smooth.dat");
		v.print();
	}

}
