package edu.uta.futureye.util;

import java.util.LinkedList;
import java.util.List;

import edu.uta.futureye.core.DOF;

/**
 * Degree Of Freedom Class
 * 自由度列表类
 * @author liuyueming
 *
 */
public class DOFList {
	protected List<DOF> dofs = new LinkedList<DOF>();

	/**
	 * @param index start from 1,2,3...
	 * @return
	 */
	public DOF at(int index) {
		if(index < 1)
			System.out.println("ERROR: DOFList index="+index);		
		return dofs.get(index-1);
	}

	public void add(DOF dof) {
		this.dofs.add(dof);
	}
	
	public int size() {
		return dofs.size();
	}
	
	public void clear() {
		dofs.clear();
	}
	
	public String toString() {
		return dofs.toString();
	}	
}
