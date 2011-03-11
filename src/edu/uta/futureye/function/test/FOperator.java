package edu.uta.futureye.function.test;

public class FOperator {
	/**
	 * coef1*f1 + coef2*f2
	 * @param coef1
	 * @param f1
	 * @param coef2
	 * @param f2
	 * @return
	 */
	public static Function LinearCombination(final double coef1, final Function f1,
			final double coef2,final Function f2) {

		Function f = new AbstractFunction(f1.getName()+" + "+f2.getName()) {
			@Override
			public void createChain() {
				chain = new PlusChain();
				PlusChain newc = null;
				if(f1.getChain() instanceof PlusChain) {
					newc = (PlusChain)f1.getChain().copy();
					newc.multiCoef(coef1);
					chain.addAllItem(newc.getAllItem());
				} else if(f1.getChain() instanceof MultiChain) {
					chain.addItem(new ItemPair(coef1,f1.getChain().copy()));
				} else {
					Exception e = new Exception("ERROR: FOperator.Plus f1");
					e.printStackTrace();
				}
				
				if(f2.getChain() instanceof PlusChain) {
					newc = (PlusChain)f2.getChain().copy();
					newc.multiCoef(coef2);
					chain.addAllItem(newc.getAllItem());
				} else if(f2.getChain() instanceof MultiChain) {
					chain.addItem(new ItemPair(coef2,f2.getChain().copy()));
				} else {
					Exception e = new Exception("ERROR: FOperator.Plus f2");
					e.printStackTrace();
				}
				chain.merge(false);
			}
		};
		return f;
	}
	
	public static Function Plus(Function f1, Function f2) {
		return LinearCombination(1.0,f1,1.0,f2) ;
	}
	
	public static Function Minus(Function f1, Function f2) {
		return LinearCombination(1.0,f1,-1.0,f2) ;
	}
	
	public static Function ScalarProduct(double s,Function f) {
		Chain c = f.getChain();
		Chain newc = null;
		if(c instanceof PlusChain) {
			newc = (PlusChain)c.copy();
			((PlusChain)newc).multiCoef(s);
		} else if(c instanceof MultiChain) {
			newc = (MultiChain)c.copy();
			newc.addItem(new ItemPair(s));
			((MultiChain)newc).merge(false);
		}
		final Chain finalc = newc;
		Function newf = new AbstractFunction(s+" * "+f.getName()) {
			@Override
			public void createChain() {
				chain = finalc;
			}
			
		};
		return newf;
		
	}

	public static Function Power(final Function f,double power) {
		MultiChain mc = new MultiChain();
		Chain c = f.getChain();
		Chain newc = null;
		if(c instanceof PlusChain) {
			newc = (PlusChain)c.copy();
			mc.addItem(new ItemPair(power,newc));
		} else if(c instanceof MultiChain) {
			mc = (MultiChain)c.copy();
			mc.multiPower(power);
		}
		final Chain finalc = mc;
		Function newf = new AbstractFunction("("+f.getName()+")^"+power) {
			@Override
			public void createChain() {
				chain = finalc;
			}
		};
		return newf;
	}
	
	public static Function Multi(Function f1,Function f2) {
		MultiChain mc = new MultiChain();
		Chain newc = null;
		Chain c1 = f1.getChain();
		if(c1 instanceof PlusChain) {
			newc = (PlusChain)c1.copy();
			mc.addItem(new ItemPair(1.0,newc));
		} else if(c1 instanceof MultiChain) {
			mc = (MultiChain)c1.copy();
		}
		
		Chain c2 = f2.getChain();
		if(c2 instanceof PlusChain) {
			newc = (PlusChain)c2.copy();
			mc.addItem(new ItemPair(1.0,newc));
		} else if(c2 instanceof MultiChain) {
			mc.addAllItem(((MultiChain)c2.copy()).getAllItem());
		} else {
			Exception e = new Exception(c2+"");
			e.printStackTrace();
		}
		
		final Chain finalc = mc;
		Function newf = new AbstractFunction("("+f1.getName()+") * ("+f2.getName()+")") {
			@Override
			public void createChain() {
				chain = finalc;
			}
		};
		return newf;
		
	}	

	public static Function Divi(Function f1,Function f2) {
		return Multi(f1,Power(f2,-1.0));
	}

	/**
	 * Ïà·´Êý
	 * @param f
	 * @return
	 */
	public static Function Opposite(final Function f) {
		Function of = new AbstractFunction() {
			@Override
			public void createChain() {
				chain = new PlusChain();
				chain.addItem(new ItemPair(-1.0, f));
				//chain.merge(true);
			}
		};
		return of;
	}
}
