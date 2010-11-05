package edu.uta.futureye.function.test;

public class FAxpb extends AbstractFunction {
	protected String varName = "x";
	protected double a;
	protected double b;

	public FAxpb(String funName,double a, double b) {
		super(funName);
		this.a = a;
		this.b = b;
	}
	
	public FAxpb(String funName,String varName, double a, double b) {
		super(funName);
		this.varName = varName;
		this.a = a;
		this.b = b;
	}
	
	@Override
	public void createChain() {
		chain = new PlusChain();
		chain.addItem(new ItemPair(a, new Variable(varName)));
		chain.addItem(new ItemPair(b));
	}
	
	public static void main(String[] args) {
		FAxpb fx = new FAxpb("F(x)=a*x+b",2.0,3.0);
		System.out.println(fx.getName());
		
		System.out.println(fx.getValue());//TODO
		
		System.out.println(fx.getValue(new Variable("x",2.0)));
		
	}

}
