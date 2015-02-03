package edu.uta.futureye.util.container;

import java.util.ArrayList;
//import scala.Function1;

public class ObjIndex extends ArrayList<Integer> {
	private static final long serialVersionUID = 6760999388784536061L;

	/**
	 * index starts from 0,1,2,3...
	 * 
	 * @param is
	 */
	public ObjIndex(Integer ...is) {
		for(Integer i : is)
			this.add(i);
	}
	
//	public <B> void foreach(Function1<Integer, B> F) {
//	    int i = 0;
//	    while (i < size()) {
//	      F.apply(get(i));
//	      i += 1;
//	    }
//	}
}
