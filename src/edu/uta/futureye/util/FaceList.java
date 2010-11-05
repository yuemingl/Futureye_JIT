package edu.uta.futureye.util;

import java.util.LinkedList;
import java.util.List;

import edu.uta.futureye.core.Face;

public class FaceList {

	protected List<Face> faces = new LinkedList<Face>();

	/**
	 * @param index start from 1,2,3...
	 * @return
	 */
	public Face at(int index) {
		if(index < 1)
			System.out.println("ERROR: FaceList index="+index);		
		return faces.get(index-1);
	}

	public void add(Face edge) {
		this.faces.add(edge);
	}
	
	public int size() {
		return faces.size();
	}
	
	public void clear() {
		faces.clear();
	}
	
	public String toString() {
		return faces.toString();
	}
}
