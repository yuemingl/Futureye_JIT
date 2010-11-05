package edu.uta.futureye.function.test;

public class ItemPair {
	public double coef;
	public Item item;
	
	public ItemPair() {
	}
	
	public ItemPair(double coef) {
		this.coef = coef;
		this.item = new GhostItem(){};
	}

	public ItemPair(double coef,Item item) {
		this.coef = coef;
		this.item = item;
	}
	
	public String toString() {
		if(item instanceof Function || item instanceof Chain)
			return coef+" * ("+item.toString()+")";
		else if(item.getName()==null)
			return ""+coef;
		return coef+" * ("+item.toString()+")";
	}
}
