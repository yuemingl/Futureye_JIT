package edu.uta.futureye.function.test;

public class FX extends AbstractFunction {
	protected String varName = "x";
	
	/**
	 * Construct function f(x) = x
	 * @param funName
	 */
	public FX(String funName) {
		super(funName);
	}
	
	public FX(String funName,String varName) {
		super(funName);
		this.varName = varName;
	}
	
	@Override
	public void createChain() {
		chain = new PlusChain();
		chain.addItem(new ItemPair(1.0, new Variable(varName)));

	}
	
	public static void main(String[] args) {
		FX fx = new FX("f(x)=x");
		System.out.println(fx.getName());
		System.out.println(fx);
		System.out.println(fx.getValue(new Variable("x",2.0)));
		
	}

}
