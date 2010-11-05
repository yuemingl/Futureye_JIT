package edu.uta.futureye.function.test;

public class ChainTest {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		MultiChain mchain = new MultiChain();
		mchain.addItem(new ItemPair(1, new Variable("x")));
		mchain.addItem(new ItemPair(2, new Variable("y")));
		mchain.addItem(new ItemPair(3, new Variable("y")));
		mchain.addItem(new ItemPair(4, new Variable("z")));
		mchain.addItem(new ItemPair(5, new Variable("z")));
		mchain.addItem(new ItemPair(6, new Variable("z")));
		mchain.addItem(new ItemPair(7, new Variable("x")));
		mchain.addItem(new ItemPair(8, new Variable("r")));
		mchain.merge(false);
		System.out.println(mchain);
	}

}
