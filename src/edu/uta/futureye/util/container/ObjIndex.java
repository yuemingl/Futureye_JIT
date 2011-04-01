package edu.uta.futureye.util.container;

import java.util.ArrayList;

public class ObjIndex extends ArrayList<Integer> {
	private static final long serialVersionUID = 6760999388784536061L;

	public ObjIndex(Integer ...is) {
		for(Integer i : is)
			this.add(i);
	}
}
