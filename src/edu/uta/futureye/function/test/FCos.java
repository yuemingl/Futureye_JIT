package edu.uta.futureye.function.test;

public class FCos extends AbstractFunction {
	Variable x = null;

	public FCos() {
		x = new Variable("x");
		setName("cos("+x.getName()+")");
	}

	public FCos(String varName) {
		x = new Variable(varName);
		setName("cos("+x.getName()+")");
	}
	
	public FCos(Variable var) {
		x = var;
		setName("cos("+x.getName()+")");
	}
	
	@Override
	public void createChain() {
		chain = new PlusChain();
		chain.addItem(new ItemPair(1.0, new ICos(x)));
	}
}
