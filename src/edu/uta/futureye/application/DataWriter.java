package edu.uta.futureye.application;

import edu.uta.futureye.algebra.SpaceVector;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.io.MeshWriter;

public class DataWriter {
	public static void writeVector(String fileName, Vector v) {
		MeshWriter w = new MeshWriter(null);
		Vector vIdx = v.copy();
		for(int i=1;i<=v.getDim();i++)
			vIdx.set(i, i);
		w.writeTechplotLine(fileName, vIdx, v);
	}
	
	public static void writeVector(String fileName, Vector index,Vector v) {
		MeshWriter w = new MeshWriter(null);
		w.writeTechplotLine(fileName, index, v);
	}
	
	public static void writeArray(String fileName, double[] ary) {
		SpaceVector v = new SpaceVector(ary);
		writeVector(fileName,v);
	}
	
	public static void writeArray(String fileName, double[] index, double[] ary) {
		SpaceVector i = new SpaceVector(index);
		SpaceVector v = new SpaceVector(ary);
		writeVector(fileName,i,v);
	}
}
