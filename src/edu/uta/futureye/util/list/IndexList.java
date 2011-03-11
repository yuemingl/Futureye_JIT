package edu.uta.futureye.util.list;

import java.util.ArrayList;

public class IndexList extends ArrayList<Integer> {
	private static final long serialVersionUID = 6760999388784536061L;

	public IndexList(Integer ...is) {
		for(Integer i : is)
			this.add(i);
	}
}
