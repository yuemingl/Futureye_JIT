package edu.uta.futureye.function.test;

import java.util.List;

public abstract class AbstractFunction extends AbstractItem implements Function {
	protected Chain chain = null;
	
	/**
	 * 匿名函数的构造函数调用父类的构造函数???
	 * if(chain == null)
			createChain();
	 */
	public AbstractFunction() {
		_name = super.toString();
	}
	
	public AbstractFunction(String name) {
		_name = name;
	}
	
	private boolean substitute(Item item, String name, Item candidate) {
		boolean rlt = false;
		if(item instanceof Chain) {
			Chain chain = (Chain)item;
			for(int i=0;i<chain.length();i++) {
				ItemPair pair = chain.getItem(i);
				if(name.equals(pair.item.getName())) {
					pair.item = candidate;
					rlt = true;
				} else {
					if(substitute(pair.item,name,candidate))
						rlt = true;
				}
			}
		} else {
			List<Item> sub = item.getSubItems();
			for(int i=0;i<sub.size();i++) {
				if(name.equals(sub.get(i).getName())) {
					sub.set(i, candidate);
					rlt = true;
				} else {
					if(substitute(sub.get(i),name,candidate)) 
						rlt = true;
				}
			}
		}
		return rlt;
	}
	
	public Function compose(final ComposePair ...pairs) {
		if(chain == null)
			createChain();
		
		final Chain chainOuter = chain;
		Function rlt = new AbstractFunction(_name) {
			@Override
			public void createChain() {
				Chain newChain = (Chain) chainOuter.copy();
				boolean found = false;
				for(int i=0;i<pairs.length;i++) {
					Function newf = (Function)pairs[i].f.copy();
					newf.setName(pairs[i].varName);
					if(substitute(newChain,pairs[i].varName,newf))
						found = true;
				}
				if(!found) {
					Exception e = new Exception("Can not substitude any items while composing function");
					e.printStackTrace();
				}
				chain = newChain; 
			}
		};

		return rlt;
	}
	
	@Override
	public double getValue(Variable... v) {
		return getValue((Item[])v);
	}
	
	@Override
	public double getValue(Item ...items) {
		if(chain == null)
			createChain();
		return chain.getValue(items);
	}
	
	@Override
	public Chain getChain() {
		if(chain == null)
			createChain();
		return chain;
	}
	
	@Override
	public void setChain(Chain chain) {
		this.chain = chain;
	}
	
	@Override
	public Item _d(String name) {
		if(chain == null)
			createChain();
		
		if(name.equals(this.getName())) {
			return new GhostItem();
		} else {
			final Chain c = (Chain) chain._d(name);
			if(c == null)
				return null; //???
			Function f = new AbstractFunction("("+_name+")_d("+name+")") {
				@Override
				public void createChain() {
					chain = c;
				}
			};
			return f;
		}
	}
	
	@Override
	public Item copy() {
		if(chain == null)
			createChain();

		final Chain c = (Chain) chain;
		Function f = new AbstractFunction(_name) {
			@Override
			public void createChain() {
				chain = (Chain) c.copy();
			}
		};
		return f;
	}
	
	@Override
	public Function expand() {
		final Chain c = COperator.ExpandChain(this.getChain());
		return new AbstractFunction(){
			@Override
			public void createChain() {
				this.chain = c;
			}
		};
	}
	
	@Override
	public int symCompairTo(Item item) {
		if(chain == null)
			createChain();		
		return this.chain.symCompairTo(item);
	}
	
	public String toString() {
		if(chain == null)
			createChain();
		return chain.toString();
	}
	
	
}