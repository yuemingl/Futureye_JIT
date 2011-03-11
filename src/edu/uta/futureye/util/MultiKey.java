package edu.uta.futureye.util;

/**
 *
 * @author liuyueming
 *
 */
public class MultiKey {
	Integer k1;
	Integer k2;
	boolean hasOrder = false;
	
	public MultiKey(int k1,int k2) {
		this.k1 = k1;
		this.k2 = k2;
	}
	
	public MultiKey(boolean hasOrder,int k1,int k2) {
		this.k1 = k1;
		this.k2 = k2;
		this.hasOrder = hasOrder;
	}
	
	
	/*****************************************************
	 * ±È½ÏË³Ðò
	 * if(obj1.hashCode() == obj2.hashCode()) {
	 *  return obj1.equals(obj2);
	 * } else {
	 * 	return false;
	 * }
	 */
	@Override
    public boolean equals(Object obj) {
    	MultiKey k = (MultiKey)obj;
    	if(this.hasOrder)
    		return (this.k1 == k.k1 && this.k2 == k.k2);
    	else {
    		return (this.k1 == k.k1 && this.k2 == k.k2) ||
    		(this.k1 == k.k2 && this.k2 == k.k1);
    	}
    }
    
	@Override
    public int hashCode() {
		//TODO
		if(this.hasOrder)
			return 70000000*k1.hashCode() + k2.hashCode();
		else
			return 70000000*(k1.hashCode() + k2.hashCode());
	}
	/*****************************************************
	 */
	
}
