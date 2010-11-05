package edu.uta.futureye.function.test;

public class GhostItem extends AbstractItem{

	@Override
	public Item copy() {
		GhostItem g = new GhostItem();
		g.setName(null);
		return g;
	}

}
