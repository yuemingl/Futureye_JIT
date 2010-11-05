package edu.uta.futureye.function.test;

public class FSin extends AbstractFunction {
	Variable x = null;

	public FSin() {
		x = new Variable("x");
		setName("sin("+x.getName()+")");
	}

	public FSin(String varName) {
		x = new Variable(varName);
		setName("sin("+x.getName()+")");
	}
	
	public FSin(Variable var) {
		x = var;
		setName("sin("+x.getName()+")");
	}
	
	@Override
	public void createChain() {
		chain = new PlusChain();
		chain.addItem(new ItemPair(1.0, new ISin(x)));
	}
}
